/* ============================================================
   Stage 1: sigma_GB(theta) for CSL boundaries, 2D LJ triangular

   Adaptive-box revision + same-box bulk reference diagnostic:
     - Box: Lx = nx*a*sqrt(Sigma), Ly = ny*a*sqrt(3*Sigma)
     - nx chosen so that Lx/2 >= max(SAFETY_D_PER_SIDE*D, MIN_LX_HALF_A*a)
       where D = a/theta_rad is the lattice-dislocation spacing
     - ny chosen so that Ly   >= MIN_SPACINGS_Y * D
     - nx AND ny forced even
     - Overlap deletion at 0.8*a0
     - relax_Lx disabled; e_ref = g_e0 from zero-pressure calibration
     - Serial compute_forces; parallelism only at sweep level

   DIAGNOSTIC (added): for every CSL entry, a perfect crystal is
   generated in the SAME box (same Lx, Ly) and its bulk energy is
   compared to g_e0. If the difference scales with Lx, then the
   use of g_e0 (computed in a 40x40 box) as reference is biasing
   the GB energy, and the fix is to use the same-box value.
   ============================================================ */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <omp.h>

#define CUTOFF          (3.0)
#define CUTOFF2         (CUTOFF * CUTOFF)
#define PI              3.141592653589793

#define OVERLAP_FRAC    0.8
#define USE_LX_RELAX    0

/* ---------- adaptive-box parameters ---------- */
#define SAFETY_D_PER_SIDE   8.0
#define MIN_LX_HALF_A       30.0
#define MIN_SPACINGS_Y      3.0

#define COARSE_NY       24
#define COARSE_MAX      800
#define COARSE_DSIG     0.2
#define COARSE_RELAX_STEPS  400
#define FINE_NY         12

#define NSX             1
#define SX_RANGE_A      0.25

#define MAX_MIN_STEPS   15000
#define FIRE_DT0        1e-3
#define FIRE_DTMAX      0.05
#define FIRE_DTMIN      1e-8
#define FIRE_DISP_CAP   0.05
#define FIRE_ALPHA0     0.1
#define FIRE_NMIN       5
#define FIRE_FTOL       1e-4
#define SHORT_FTOL      1e-3

#define LX_RELAX_STEPS  8
#define LX_RELAX_DELTA  0.002

#define MAX_MN          25
#define MAX_SIGMA       450
#define MIN_ANGLE_DEG   4.5
#define MAX_ANGLE_DEG   40.0
#define MIN_MULT        2

static double g_a0   = 0.0;
static double g_e0   = 0.0;
static double g_rho0 = 0.0;

static double SF_U_C = 0.0, SF_F_C = 0.0;
static void init_shifted_force(void) {
    double rc = CUTOFF;
    double r2i = 1.0/(rc*rc), r6i = r2i*r2i*r2i, r12i = r6i*r6i;
    SF_U_C = 4.0*(r12i - r6i);
    SF_F_C = 48.0/rc*(r12i - 0.5*r6i);
}

typedef struct {
    int N;
    double *x, *y, *vx, *vy, *fx, *fy;
    double Lx, Ly, origin_x, origin_y;
} Atoms;

static void *xmalloc(size_t s){void*p=malloc(s);if(!p){fprintf(stderr,"OOM\n");exit(1);}return p;}
static void *xcalloc(size_t n,size_t s){void*p=calloc(n,s);if(!p){fprintf(stderr,"OOM\n");exit(1);}return p;}
static void free_atoms(Atoms*a){
    free(a->x);free(a->y);free(a->vx);free(a->vy);free(a->fx);free(a->fy);
}

static _Thread_local int *tls_head      = NULL;
static _Thread_local int *tls_lscl      = NULL;
static _Thread_local int  tls_ncell_cap = 0;
static _Thread_local int  tls_N_cap     = 0;

static void tls_ensure(int ncell, int N) {
    if (ncell > tls_ncell_cap) {
        free(tls_head);
        tls_head = xmalloc((size_t)ncell * sizeof(int));
        tls_ncell_cap = ncell;
    }
    if (N > tls_N_cap) {
        free(tls_lscl);
        tls_lscl = xmalloc((size_t)N * sizeof(int));
        tls_N_cap = N;
    }
}

static inline void min_image(const Atoms *a, double xi,double yi,double xj,double yj,
                             double *dx,double *dy){
    double ddx=xi-xj, ddy=yi-yj;
    ddx -= a->Lx*round(ddx/a->Lx);
    ddy -= a->Ly*round(ddy/a->Ly);
    *dx=ddx; *dy=ddy;
}
static inline void wrap_pbc(const Atoms *a, double *x,double *y){
    double sx=(*x-a->origin_x)/a->Lx, sy=(*y-a->origin_y)/a->Ly;
    sx-=floor(sx); sy-=floor(sy);
    *x=a->origin_x+sx*a->Lx; *y=a->origin_y+sy*a->Ly;
}
static inline void lj_force_pe(double dx,double dy,
                               double *fx,double *fy,double *pe){
    double r2=dx*dx+dy*dy;
    if(r2>=CUTOFF2||r2<1e-4){*fx=0;*fy=0;*pe=0;return;}
    double r=sqrt(r2), r2i=1.0/r2, r6i=r2i*r2i*r2i, r12i=r6i*r6i;
    double U_raw = 4.0*(r12i - r6i);
    double F_raw = 48.0/r*(r12i - 0.5*r6i);
    double U_sf  = U_raw - SF_U_C + SF_F_C*(r - CUTOFF);
    double F_sf  = F_raw - SF_F_C;
    double f_or  = F_sf / r;
    *fx = f_or*dx; *fy = f_or*dy; *pe = U_sf;
}

static void compute_forces_allpairs(Atoms *a, double *pe_out){
    int N=a->N;
    double *fx=a->fx, *fy=a->fy;
    for(int i=0;i<N;i++){fx[i]=0;fy[i]=0;}
    double pe=0.0;
    for(int i=0;i<N;i++){
        for(int j=i+1;j<N;j++){
            double dx,dy;
            min_image(a,a->x[i],a->y[i],a->x[j],a->y[j],&dx,&dy);
            double fxp,fyp,pep;
            lj_force_pe(dx,dy,&fxp,&fyp,&pep);
            fx[i]+=fxp; fy[i]+=fyp;
            fx[j]-=fxp; fy[j]-=fyp;
            pe += pep;
        }
    }
    *pe_out=pe;
}

void compute_forces(Atoms *a, double *pe_out){
    int N=a->N;
    double *x=a->x,*y=a->y,*fx=a->fx,*fy=a->fy;
    double Lx=a->Lx, Ly=a->Ly, ox=a->origin_x, oy=a->origin_y;

    for(int i=0;i<N;i++){fx[i]=0;fy[i]=0;}
    if (N == 0) { *pe_out = 0.0; return; }

    int nx=(int)(Lx/CUTOFF); if(nx<1)nx=1;
    int ny=(int)(Ly/CUTOFF); if(ny<1)ny=1;

    if (nx < 3 || ny < 3) {
        compute_forces_allpairs(a, pe_out);
        return;
    }

    double cxs=Lx/nx, cys=Ly/ny;
    int ncell=nx*ny;

    tls_ensure(ncell, N);
    int *head = tls_head;
    int *lscl = tls_lscl;

    for(int c=0;c<ncell;c++) head[c]=-1;
    for(int i=0;i<N;i++){
        int cx=(int)((x[i]-ox)/cxs);
        int cy=(int)((y[i]-oy)/cys);
        cx=((cx%nx)+nx)%nx; cy=((cy%ny)+ny)%ny;
        int cell=cy*nx+cx;
        lscl[i]=head[cell]; head[cell]=i;
    }

    double pe=0.0;
    for(int i=0;i<N;i++){
        int cx=(int)((x[i]-ox)/cxs);
        int cy=(int)((y[i]-oy)/cys);
        cx=((cx%nx)+nx)%nx; cy=((cy%ny)+ny)%ny;
        for(int di=-1;di<=1;di++){
            int iix=((cx+di)%nx+nx)%nx;
            for(int dj=-1;dj<=1;dj++){
                int jjy=((cy+dj)%ny+ny)%ny;
                int cell=jjy*nx+iix;
                for(int j=head[cell];j!=-1;j=lscl[j]){
                    if(j<=i) continue;
                    double dx,dy;
                    min_image(a,x[i],y[i],x[j],y[j],&dx,&dy);
                    double fxp,fyp,pep;
                    lj_force_pe(dx,dy,&fxp,&fyp,&pep);
                    fx[i]+=fxp; fy[i]+=fyp;
                    fx[j]-=fxp; fy[j]-=fyp;
                    pe += pep;
                }
            }
        }
    }
    *pe_out=pe;
}

double compute_total_pe(Atoms *a){
    int N=a->N;
    double Lx=a->Lx, Ly=a->Ly, ox=a->origin_x, oy=a->origin_y;

    int nx=(int)(Lx/CUTOFF); if(nx<1)nx=1;
    int ny=(int)(Ly/CUTOFF); if(ny<1)ny=1;

    if (nx < 3 || ny < 3) {
        double pe=0.0;
        for(int i=0;i<N;i++){
            for(int j=i+1;j<N;j++){
                double dx,dy;
                min_image(a,a->x[i],a->y[i],a->x[j],a->y[j],&dx,&dy);
                double fxp,fyp,pep;
                lj_force_pe(dx,dy,&fxp,&fyp,&pep);
                pe += pep;
            }
        }
        return pe;
    }

    double cxs=Lx/nx, cys=Ly/ny;
    int ncell=nx*ny;

    tls_ensure(ncell, N);
    int *head = tls_head;
    int *lscl = tls_lscl;

    for(int c=0;c<ncell;c++) head[c]=-1;
    for(int i=0;i<N;i++){
        int cx=(int)((a->x[i]-ox)/cxs);
        int cy=(int)((a->y[i]-oy)/cys);
        cx=((cx%nx)+nx)%nx; cy=((cy%ny)+ny)%ny;
        int cell=cy*nx+cx;
        lscl[i]=head[cell]; head[cell]=i;
    }
    double pe=0.0;
    for(int i=0;i<N;i++){
        int cx=(int)((a->x[i]-ox)/cxs);
        int cy=(int)((a->y[i]-oy)/cys);
        cx=((cx%nx)+nx)%nx; cy=((cy%ny)+ny)%ny;
        for(int di=-1;di<=1;di++){
            int iix=((cx+di)%nx+nx)%nx;
            for(int dj=-1;dj<=1;dj++){
                int jjy=((cy+dj)%ny+ny)%ny;
                int cell=jjy*nx+iix;
                for(int j=head[cell];j!=-1;j=lscl[j]){
                    if(j<=i) continue;
                    double dx,dy;
                    min_image(a,a->x[i],a->y[i],a->x[j],a->y[j],&dx,&dy);
                    double fxp,fyp,pep;
                    lj_force_pe(dx,dy,&fxp,&fyp,&pep);
                    pe += pep;
                }
            }
        }
    }
    return pe;
}

static void remove_overlaps(Atoms *A, double rcut){
    int N=A->N;
    if(N<2) return;
    double Lx=A->Lx, Ly=A->Ly, ox=A->origin_x, oy=A->origin_y;
    double rmin2=rcut*rcut;

    int nx=(int)(Lx/rcut); if(nx<1)nx=1;
    int ny=(int)(Ly/rcut); if(ny<1)ny=1;
    double cxs=Lx/nx, cys=Ly/ny;
    int ncell=nx*ny;

    int *head=xmalloc(ncell*sizeof(int));
    int *lscl=xmalloc(N*sizeof(int));
    for(int c=0;c<ncell;c++) head[c]=-1;
    for(int i=0;i<N;i++){
        int cx=(int)((A->x[i]-ox)/cxs);
        int cy=(int)((A->y[i]-oy)/cys);
        cx=((cx%nx)+nx)%nx; cy=((cy%ny)+ny)%ny;
        int cell=cy*nx+cx;
        lscl[i]=head[cell]; head[cell]=i;
    }

    bool *keep=xmalloc(N*sizeof(bool));
    for(int i=0;i<N;i++) keep[i]=true;

    for(int i=0;i<N;i++){
        if(!keep[i]) continue;
        int cx=(int)((A->x[i]-ox)/cxs);
        int cy=(int)((A->y[i]-oy)/cys);
        cx=((cx%nx)+nx)%nx; cy=((cy%ny)+ny)%ny;
        for(int di=-1;di<=1;di++){
            int ii=((cx+di)%nx+nx)%nx;
            for(int dj=-1;dj<=1;dj++){
                int jj=((cy+dj)%ny+ny)%ny;
                int cell=jj*nx+ii;
                for(int j=head[cell];j!=-1;j=lscl[j]){
                    if(j<=i) continue;
                    if(!keep[j]) continue;
                    double dx,dy;
                    min_image(A,A->x[i],A->y[i],A->x[j],A->y[j],&dx,&dy);
                    if(dx*dx+dy*dy < rmin2) {
                        if ( ((i + j) & 1) == 0 ) keep[j]=false;
                        else                       keep[i]=false;
                    }
                }
            }
        }
    }

    int Nn=0;
    for(int i=0;i<N;i++) if(keep[i]){
        A->x[Nn]=A->x[i]; A->y[Nn]=A->y[i]; Nn++;
    }
    A->N=Nn;
    free(head); free(lscl); free(keep);
}

static double min_pair_distance(Atoms *a){
    int N=a->N;
    if(N<2) return 1e30;
    double rmin2=1e30;
    for(int i=0;i<N;i++){
        for(int j=i+1;j<N;j++){
            double dx,dy;
            min_image(a,a->x[i],a->y[i],a->x[j],a->y[j],&dx,&dy);
            double r2=dx*dx+dy*dy;
            if(r2<rmin2) rmin2=r2;
        }
    }
    return sqrt(rmin2);
}

static int fire_minimise_impl(Atoms *a, int max_steps, double ftol,
                              int quiet, double *pe_out){
    (void)quiet;
    int N=a->N;
    double dt=FIRE_DT0, alpha=FIRE_ALPHA0;
    int n_pos=0;
    double pe;
    compute_forces(a,&pe);
    for(int i=0;i<N;i++){a->vx[i]=0;a->vy[i]=0;}

    int converged=0;
    for(int step=0;step<max_steps;step++){
        for(int i=0;i<N;i++){
            a->vx[i] += 0.5*dt*a->fx[i];
            a->vy[i] += 0.5*dt*a->fy[i];
        }
        for(int i=0;i<N;i++){
            double dx = dt*a->vx[i];
            double dy = dt*a->vy[i];
            double d2 = dx*dx + dy*dy;
            double cap2 = FIRE_DISP_CAP*FIRE_DISP_CAP;
            if (d2 > cap2) {
                double scale = sqrt(cap2/d2);
                dx *= scale; dy *= scale;
                a->vx[i] *= scale; a->vy[i] *= scale;
            }
            a->x[i] += dx;
            a->y[i] += dy;
        }
        for(int i=0;i<N;i++) wrap_pbc(a,&a->x[i],&a->y[i]);
        compute_forces(a,&pe);
        for(int i=0;i<N;i++){
            a->vx[i] += 0.5*dt*a->fx[i];
            a->vy[i] += 0.5*dt*a->fy[i];
        }
        double P=0, v2=0, f2=0;
        for(int i=0;i<N;i++){
            P  += a->fx[i]*a->vx[i] + a->fy[i]*a->vy[i];
            v2 += a->vx[i]*a->vx[i] + a->vy[i]*a->vy[i];
            f2 += a->fx[i]*a->fx[i] + a->fy[i]*a->fy[i];
        }
        double vn=sqrt(v2), fn=sqrt(f2);
        if(P>0.0){
            n_pos++;
            if(n_pos>FIRE_NMIN){
                dt = dt*1.1;
                if(dt>FIRE_DTMAX) dt=FIRE_DTMAX;
                alpha *= 0.99;
                if(alpha<0.01) alpha=0.01;
            }
            if(vn>1e-14 && fn>1e-14){
                double scale = alpha*vn/fn;
                for(int i=0;i<N;i++){
                    a->vx[i] = (1.0-alpha)*a->vx[i] + scale*a->fx[i];
                    a->vy[i] = (1.0-alpha)*a->vy[i] + scale*a->fy[i];
                }
            }
        } else {
            n_pos = 0;
            dt *= 0.5;
            if(dt < FIRE_DTMIN) dt = FIRE_DTMIN;
            alpha = FIRE_ALPHA0;
            for(int i=0;i<N;i++){a->vx[i]=0;a->vy[i]=0;}
        }
        if((step & 127) == 0){
            double fmax2=0;
            for(int i=0;i<N;i++){
                double ff = a->fx[i]*a->fx[i] + a->fy[i]*a->fy[i];
                if(ff>fmax2) fmax2=ff;
            }
            if(sqrt(fmax2) < ftol){ converged=1; break; }
        }
    }
    if(!converged){
        double fmax2=0;
        for(int i=0;i<N;i++){
            double ff = a->fx[i]*a->fx[i] + a->fy[i]*a->fy[i];
            if(ff>fmax2) fmax2=ff;
        }
        if(sqrt(fmax2) < ftol) converged=1;
    }
    if(pe_out) *pe_out = pe;
    return converged;
}

int fire_minimise(Atoms *a, int max_steps, double *pe_out){
    return fire_minimise_impl(a, max_steps, FIRE_FTOL, 0, pe_out);
}
static int fire_minimise_short(Atoms *a, int max_steps, double *pe_out){
    return fire_minimise_impl(a, max_steps, SHORT_FTOL, 1, pe_out);
}

static void relax_Lx(Atoms *a, double *pe_io){
    double pe_cur = *pe_io;
    double Lx_save = a->Lx;
    double ox_save = a->origin_x;
    double *x_save = xmalloc(a->N * sizeof(double));
    memcpy(x_save, a->x, a->N * sizeof(double));

    double scale = 1.0 + LX_RELAX_DELTA;
    for(int i=0;i<a->N;i++) a->x[i] *= scale;
    a->Lx *= scale; a->origin_x *= scale;
    for(int i=0;i<a->N;i++) wrap_pbc(a,&a->x[i],&a->y[i]);
    double pe_new = compute_total_pe(a);

    if (pe_new < pe_cur) { *pe_io = pe_new; free(x_save); return; }

    a->Lx = Lx_save; a->origin_x = ox_save;
    memcpy(a->x, x_save, a->N * sizeof(double));
    scale = 1.0 - LX_RELAX_DELTA;
    for(int i=0;i<a->N;i++) a->x[i] *= scale;
    a->Lx *= scale; a->origin_x *= scale;
    for(int i=0;i<a->N;i++) wrap_pbc(a,&a->x[i],&a->y[i]);
    pe_new = compute_total_pe(a);

    if (pe_new < pe_cur) { *pe_io = pe_new; free(x_save); return; }

    a->Lx = Lx_save; a->origin_x = ox_save;
    memcpy(a->x, x_save, a->N * sizeof(double));
    for(int i=0;i<a->N;i++) wrap_pbc(a,&a->x[i],&a->y[i]);
    free(x_save);
}

typedef struct {
    int m, n;
    int Sigma;
    int Sigma_true;
    double theta_deg;
    double Lx_a, Ly_a;
    int nx, ny;
} CSL_Entry;

static int gcd_int(int a,int b){
    while(b){int t=a%b;a=b;b=t;}
    return a;
}

static void choose_multipliers(int Sig, double theta_rad, int *nx_out, int *ny_out){
    double D_over_a = 1.0 / theta_rad;

    double Lx_half_min_a = SAFETY_D_PER_SIDE * D_over_a;
    if(Lx_half_min_a < MIN_LX_HALF_A) Lx_half_min_a = MIN_LX_HALF_A;
    double Lx_min_a = 2.0 * Lx_half_min_a;
    double Lx_base_a = sqrt((double)Sig);
    int nx = (int)ceil(Lx_min_a / Lx_base_a);
    if(nx < MIN_MULT) nx = MIN_MULT;
    if(nx % 2 != 0) nx++;

    double Ly_min_a = MIN_SPACINGS_Y * D_over_a;
    double Ly_base_a = sqrt(3.0*(double)Sig);
    int ny = (int)ceil(Ly_min_a / Ly_base_a);
    if(ny < MIN_MULT) ny = MIN_MULT;
    if(ny % 2 != 0) ny++;

    *nx_out = nx; *ny_out = ny;
}

int generate_csl(CSL_Entry *tab,int cap){
    int count=0;
    for(int m=1;m<=MAX_MN;m++){
        for(int n=0;n<m;n++){
            if(n==0 && m==1) continue;
            if(gcd_int(m,n)!=1) continue;
            int Sig = m*m + m*n + n*n;
            if(Sig > MAX_SIGMA) continue;

            double theta_deg = 2.0*atan(sqrt(3.0)*n / (2.0*m + n)) * 180.0/PI;
            if(theta_deg < MIN_ANGLE_DEG || theta_deg > MAX_ANGLE_DEG) continue;
            double theta_rad = theta_deg * PI / 180.0;

            int Sig_true = (Sig % 3 == 0) ? (Sig/3) : Sig;

            int dup=-1;
            for(int i=0;i<count;i++){
                if(fabs(tab[i].theta_deg - theta_deg) < 0.05){dup=i;break;}
            }
            if(dup>=0){
                if(Sig < tab[dup].Sigma){
                    int nx, ny;
                    choose_multipliers(Sig, theta_rad, &nx, &ny);
                    tab[dup].m=m; tab[dup].n=n;
                    tab[dup].Sigma=Sig;
                    tab[dup].Sigma_true=Sig_true;
                    tab[dup].theta_deg=theta_deg;
                    tab[dup].Lx_a = sqrt((double)Sig);
                    tab[dup].Ly_a = sqrt(3.0*(double)Sig);
                    tab[dup].nx = nx; tab[dup].ny = ny;
                }
                continue;
            }

            double Lx_a = sqrt((double)Sig);
            double Ly_a = sqrt(3.0*(double)Sig);
            int nx, ny;
            choose_multipliers(Sig, theta_rad, &nx, &ny);

            tab[count].m=m; tab[count].n=n;
            tab[count].Sigma=Sig;
            tab[count].Sigma_true=Sig_true;
            tab[count].theta_deg=theta_deg;
            tab[count].Lx_a=Lx_a;
            tab[count].Ly_a=Ly_a;
            tab[count].nx=nx; tab[count].ny=ny;
            count++;
            if(count>=cap) return count;
        }
    }

    for(int i=0;i<count-1;i++){
        for(int j=i+1;j<count;j++){
            if(tab[j].theta_deg < tab[i].theta_deg){
                CSL_Entry t = tab[i];
                tab[i] = tab[j];
                tab[j] = t;
            }
        }
    }

    for(int i=0;i<count-1;i++){
        if(fabs(tab[i].theta_deg - tab[i+1].theta_deg) < 1e-6){
            fprintf(stderr,"DUPLICATE: theta=%.6f at indices %d,%d\n",
                    tab[i].theta_deg, i, i+1);
        }
    }

    return count;
}

Atoms generate_bicrystal(const CSL_Entry *csl, double shift_x_a, double shift_y_a){
    double a = g_a0;
    double theta = csl->theta_deg * PI/180.0;
    double ch = cos(0.5*theta), sh = sin(0.5*theta);

    double Lx = (double)csl->nx * csl->Lx_a * a;
    double Ly = (double)csl->ny * csl->Ly_a * a;

    Atoms A;
    A.N=0;
    A.Lx=Lx; A.Ly=Ly;
    A.origin_x=-0.5*Lx; A.origin_y=-0.5*Ly;

    double sx = 0.5 * shift_x_a * a;
    double sy = shift_y_a * a;

    int max_atoms = (int)(3.0*Lx*Ly/(a*a*sqrt(3.0))) + 1000;
    A.x=xmalloc(max_atoms*sizeof(double));
    A.y=xmalloc(max_atoms*sizeof(double));

    int oversample = 3;
    int Ngen = (int)(Lx/a) * oversample + 8;
    int Mgen = (int)(Ly/(a*sqrt(3.0)/2.0)) * oversample + 8;
    double dy_row = a*sqrt(3.0)/2.0;
    double eps_x = 1e-9 * Lx;
    double eps_y = 1e-9 * Ly;

    for(int r=-Mgen;r<=Mgen;r++){
        double y0 = (r + 0.5)*dy_row;
        double xoff = ((r%2)+2)%2 * 0.5*a;
        for(int c=-Ngen;c<=Ngen;c++){
            double x0 = c*a + xoff;
            double xr = ch*x0 - sh*y0 - sx;
            double yr = sh*x0 + ch*y0 + sy;
            if(xr < -eps_x && xr > -Lx/2.0 - eps_x &&
               yr >= -Ly/2.0 - eps_y && yr < Ly/2.0 - eps_y){
                wrap_pbc(&A,&xr,&yr);
                if(A.N >= max_atoms){fprintf(stderr,"ERR max_atoms 1\n");exit(1);}
                A.x[A.N]=xr; A.y[A.N]=yr; A.N++;
            }
        }
    }
    for(int r=-Mgen;r<=Mgen;r++){
        double y0 = (r + 0.5)*dy_row;
        double xoff = ((r%2)+2)%2 * 0.5*a;
        for(int c=-Ngen;c<=Ngen;c++){
            double x0 = c*a + xoff;
            double xr = ch*x0 + sh*y0 + sx;
            double yr = -sh*x0 + ch*y0;
            if(xr >= -eps_x && xr < Lx/2.0 - eps_x &&
               yr >= -Ly/2.0 - eps_y && yr < Ly/2.0 - eps_y){
                wrap_pbc(&A,&xr,&yr);
                if(A.N >= max_atoms){fprintf(stderr,"ERR max_atoms 2\n");exit(1);}
                A.x[A.N]=xr; A.y[A.N]=yr; A.N++;
            }
        }
    }

    remove_overlaps(&A, OVERLAP_FRAC * g_a0);

    A.vx=xcalloc(A.N,sizeof(double));
    A.vy=xcalloc(A.N,sizeof(double));
    A.fx=xmalloc(A.N*sizeof(double));
    A.fy=xmalloc(A.N*sizeof(double));
    return A;
}

Atoms generate_perfect_at(int Lx_atoms, int Ly_rows, double a){
    double dy_row = a*sqrt(3.0)/2.0;
    double Lx = Lx_atoms * a;
    double Ly = Ly_rows * dy_row;

    Atoms A;
    A.N=0;
    A.Lx=Lx; A.Ly=Ly;
    A.origin_x=-0.5*Lx; A.origin_y=-0.5*Ly;

    A.x=xmalloc((size_t)Lx_atoms*Ly_rows*sizeof(double)+64);
    A.y=xmalloc((size_t)Lx_atoms*Ly_rows*sizeof(double)+64);
    for(int r=0;r<Ly_rows;r++){
        double y_pos = A.origin_y + (r + 0.5)*dy_row;
        double xoff = (r%2) * 0.5*a;
        for(int c=0;c<Lx_atoms;c++){
            double x_pos = A.origin_x + c*a + xoff;
            wrap_pbc(&A,&x_pos,&y_pos);
            A.x[A.N]=x_pos; A.y[A.N]=y_pos; A.N++;
        }
    }
    A.vx=xcalloc(A.N,sizeof(double));
    A.vy=xcalloc(A.N,sizeof(double));
    A.fx=xmalloc(A.N*sizeof(double));
    A.fy=xmalloc(A.N*sizeof(double));
    return A;
}

static void find_a0_and_e0(void){
    const int LX = 40, LY = 40;
    double a_lo = 1.06, a_hi = 1.18;
    int    Nscan = 61;
    double a_best = a_lo;
    double e_best = 1e30;

    for(int k=0;k<Nscan;k++){
        double a = a_lo + (a_hi - a_lo) * k / (Nscan - 1);
        Atoms P = generate_perfect_at(LX, LY, a);
        double pe = compute_total_pe(&P);
        double e = pe / P.N;
        if(e < e_best){ e_best = e; a_best = a; }
        free_atoms(&P);
    }

    double da = 0.001;
    Atoms Pa = generate_perfect_at(LX, LY, a_best - da);
    Atoms Pb = generate_perfect_at(LX, LY, a_best);
    Atoms Pc = generate_perfect_at(LX, LY, a_best + da);
    double ea = compute_total_pe(&Pa) / Pa.N;
    double eb = compute_total_pe(&Pb) / Pb.N;
    double ec = compute_total_pe(&Pc) / Pc.N;
    free_atoms(&Pa); free_atoms(&Pb); free_atoms(&Pc);

    double denom = (ea - 2.0*eb + ec);
    if (fabs(denom) > 1e-15) {
        double offset = 0.5 * (ea - ec) / denom;
        if (fabs(offset) < 1.0) {
            a_best += offset * da;
            Atoms Pr = generate_perfect_at(LX, LY, a_best);
            e_best = compute_total_pe(&Pr) / Pr.N;
            free_atoms(&Pr);
        }
    }

    g_a0   = a_best;
    g_e0   = e_best;
    g_rho0 = 2.0 / (a_best * a_best * sqrt(3.0));
}

double compute_sigma_gb(const CSL_Entry *csl, double e_ref,
                        double *best_sx_out, double *best_sy_out,
                        int *N_out, double *pe_out,
                        int *converged_out, double *dmin_out)
{
    double dy_range = csl->Ly_a;

    int n_coarse = (int)ceil(dy_range / COARSE_DSIG);
    if (n_coarse < COARSE_NY)  n_coarse = COARSE_NY;
    if (n_coarse > COARSE_MAX) n_coarse = COARSE_MAX;
    double coarse_spacing = dy_range / (double)n_coarse;
    double fine_range     = coarse_spacing;

    double *sys = xmalloc((size_t)n_coarse*sizeof(double));
    double *pes = xmalloc((size_t)n_coarse*sizeof(double));
    for(int j=0;j<n_coarse;j++)
        sys[j] = (j/(double)n_coarse - 0.5)*dy_range;

    #pragma omp parallel for schedule(dynamic)
    for(int j=0;j<n_coarse;j++){
        Atoms A = generate_bicrystal(csl, 0.0, sys[j]);
        double pe;
        fire_minimise_short(&A, COARSE_RELAX_STEPS, &pe);
        pes[j] = (pe - A.N*e_ref) / (2.0*A.Ly);
        free_atoms(&A);
    }
    int best=0;
    for(int j=1;j<n_coarse;j++) if(pes[j]<pes[best]) best=j;
    double best_sy_A = sys[best];
    double best_e    = pes[best];
    free(sys); free(pes);

    double *fys = xmalloc(FINE_NY*sizeof(double));
    for(int j=0;j<FINE_NY;j++)
        fys[j] = best_sy_A + ((j/(double)(FINE_NY-1)) - 0.5)*2.0*fine_range;

    double *fpes = xmalloc(FINE_NY*sizeof(double));
    #pragma omp parallel for schedule(dynamic)
    for(int j=0;j<FINE_NY;j++){
        Atoms A = generate_bicrystal(csl, 0.0, fys[j]);
        double pe;
        fire_minimise_short(&A, COARSE_RELAX_STEPS, &pe);
        fpes[j] = (pe - A.N*e_ref)/(2.0*A.Ly);
        free_atoms(&A);
    }
    for(int j=0;j<FINE_NY;j++){
        if(fpes[j] < best_e){ best_e=fpes[j]; best_sy_A=fys[j]; }
    }
    free(fys); free(fpes);

    double sx_min = -SX_RANGE_A, sx_max = SX_RANGE_A;
    double best_sx = 0.0;
    double best_sy = best_sy_A;

    int n_sx = NSX;
    int n_sy = FINE_NY;
    double *sx_arr = xmalloc(n_sx * sizeof(double));
    double *sy_arr = xmalloc(n_sy * sizeof(double));
    for(int i=0;i<n_sx;i++){
        if (n_sx == 1) sx_arr[i] = 0.0;
        else           sx_arr[i] = sx_min + (sx_max - sx_min) * i / (double)(n_sx - 1);
    }
    for(int j=0;j<n_sy;j++){
        if (n_sy == 1) sy_arr[j] = best_sy_A;
        else           sy_arr[j] = best_sy_A + ((j/(double)(n_sy-1)) - 0.5)*2.0*fine_range;
    }
    double *pes2 = xmalloc((size_t)n_sx * n_sy * sizeof(double));
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for(int i=0;i<n_sx;i++){
        for(int j=0;j<n_sy;j++){
            Atoms A = generate_bicrystal(csl, sx_arr[i], sy_arr[j]);
            double pe;
            fire_minimise_short(&A, COARSE_RELAX_STEPS, &pe);
            pes2[i*n_sy+j] = (pe - A.N*e_ref)/(2.0*A.Ly);
            free_atoms(&A);
        }
    }
    for(int i=0;i<n_sx;i++){
        for(int j=0;j<n_sy;j++){
            double e = pes2[i*n_sy+j];
            if(e < best_e){
                best_e  = e;
                best_sx = sx_arr[i];
                best_sy = sy_arr[j];
            }
        }
    }
    free(sx_arr); free(sy_arr); free(pes2);

    Atoms A = generate_bicrystal(csl, best_sx, best_sy);

    if(A.N < 10){
        fprintf(stderr,"  WARN: theta=%.4f (Sig=%d) N=%d -- skipping\n",
                csl->theta_deg, csl->Sigma, A.N);
        *best_sx_out = best_sx;
        *best_sy_out = best_sy;
        *N_out = A.N;
        *pe_out = 0.0;
        *converged_out = 0;
        *dmin_out = 0.0;
        free_atoms(&A);
        return NAN;
    }

    /* ============================================================
       DIAGNOSTIC: same-box bulk reference.
       Generate a perfect crystal in the SAME box (Lx, Ly) as the
       bicrystal and compare its energy per atom to g_e0.
       If the difference scales with Lx, then the use of g_e0 is
       biasing the GB energy and the fix is to use e_bulk_same.
       ============================================================ */
    {
        int ref_cols = (int)round(A.Lx / g_a0);
        int ref_rows = (int)round(A.Ly / (g_a0 * sqrt(3.0)/2.0));
        if(ref_rows % 2 != 0) ref_rows++;
        Atoms R = generate_perfect_at(ref_cols, ref_rows, g_a0);
        double e_bulk_same = compute_total_pe(&R) / R.N;
        printf("  [REF] theta=%8.4f  Lx/a=%7.1f  Lx/(2D)=%6.2f  "
               "e_bulk_same=%.9f  g_e0=%.9f  diff=%+.2e\n",
               csl->theta_deg,
               A.Lx/g_a0,
               0.5*(A.Lx/g_a0) / (1.0 / (csl->theta_deg * PI/180.0)),
               e_bulk_same, g_e0, e_bulk_same - g_e0);
        free_atoms(&R);
    }

    double pe_final;
    fire_minimise(&A, MAX_MIN_STEPS, &pe_final);

#if USE_LX_RELAX
    for(int k=0;k<LX_RELAX_STEPS;k++){
        double pe_before = pe_final;
        relax_Lx(&A, &pe_final);
        if (fabs(pe_final - pe_before) < 1e-6) break;
        fire_minimise(&A, MAX_MIN_STEPS, &pe_final);
    }
#endif

    double sigma = (pe_final - A.N*e_ref) / (2.0*A.Ly);
    double dmin = min_pair_distance(&A);

    *best_sx_out = best_sx;
    *best_sy_out = best_sy;
    *N_out = A.N;
    *pe_out = pe_final;
    *converged_out = 1;
    *dmin_out = dmin;
    free_atoms(&A);
    return sigma;
}

int main(void){
    printf("==========================================================\n");
    printf("Stage 1: sigma_GB(theta), 2D LJ triangular (adaptive box)\n");
    printf("==========================================================\n");

#if defined(_OPENMP) && _OPENMP >= 200805
    omp_set_max_active_levels(1);
#else
    omp_set_nested(0);
#endif

    init_shifted_force();
    printf("Shifted-force LJ: U(rc)=%.6f, F(rc)=%.6f\n", SF_U_C, SF_F_C);

    printf("\n[Calibration] finding a0 at P=0 (minimize E(a))...\n");
    find_a0_and_e0();
    printf("  a0       = %.6f\n", g_a0);
    printf("  rho0     = %.6f\n", g_rho0);
    printf("  e0 (P=0) = %.6f eps/atom\n", g_e0);

    printf("\nConfiguration:\n");
    printf("  a = a0 = %.6f, cutoff=%.3f, overlap_frac=%.2f (-> %.3f sigma)\n",
           g_a0, CUTOFF, OVERLAP_FRAC, OVERLAP_FRAC * g_a0);
    printf("  FIRE: FTOL=%.0e, DT0=%.0e, DTMIN=%.0e, DISP_CAP=%.3f\n",
           FIRE_FTOL, FIRE_DT0, FIRE_DTMIN, FIRE_DISP_CAP);
    printf("  Adaptive box: Lx/2 >= max(%.1f*D, %.1f*a), Ly >= %.1f*D\n",
           SAFETY_D_PER_SIDE, MIN_LX_HALF_A, MIN_SPACINGS_Y);
    printf("  nx AND ny forced even\n");
    printf("  Lx relaxation: %s\n", USE_LX_RELAX ? "ON" : "OFF");
    printf("  OpenMP threads: %d (parallelism only at sweep level)\n",
           omp_get_max_threads());

    double e_ref = g_e0;
    printf("\n[Reference] e_ref = %.6f eps/atom\n", e_ref);

    CSL_Entry *csl = xmalloc(MAX_MN*MAX_MN*sizeof(CSL_Entry));
    int n_csl = generate_csl(csl, MAX_MN*MAX_MN);
    printf("\nGenerated %d CSL entries\n", n_csl);
    printf("  %-8s %-6s %-6s %-8s %-8s %-8s %-8s\n",
           "theta", "Sig", "Sig_t", "(m,n)", "nx,ny", "Lx/a", "Lx/(2D)");
    printf("  -----------------------------------------------------------------\n");
    for(int i=0;i<n_csl;i++){
        double theta_rad = csl[i].theta_deg * PI / 180.0;
        double D_over_a  = 1.0 / theta_rad;
        double Lx_a = csl[i].nx * csl[i].Lx_a;
        double Lx_half_over_D = (0.5*Lx_a) / D_over_a;
        printf("  %-8.4f %-6d %-6d (%d,%d)    %d,%d     %-8.2f %-8.2f\n",
               csl[i].theta_deg, csl[i].Sigma, csl[i].Sigma_true,
               csl[i].m, csl[i].n, csl[i].nx, csl[i].ny,
               Lx_a, Lx_half_over_D);
    }

    FILE *out = fopen("sigma_gb_csl.csv","w");
    if(!out){fprintf(stderr,"cannot open csv\n");return 1;}
    fprintf(out,"theta_deg,sigma_gb,Sigma,Sigma_true,m,n,nx,ny,N_atoms,"
                "best_sx_a,best_sy_a,pe_final,converged,dmin\n");

    printf("\n[Sweep]\n");
    printf("  %-8s %-6s %-6s %-12s %-10s %-10s %-8s %-8s\n",
           "theta", "Sig", "Sig_t", "sigma_gb", "sx/a", "sy/a", "N", "dmin");
    printf("  ---------------------------------------------------------------------------\n");

    for(int i=0;i<n_csl;i++){
        double sx, sy, pe, dmin;
        int N, conv;
        double sig = compute_sigma_gb(&csl[i], e_ref, &sx, &sy, &N, &pe, &conv, &dmin);

        if(isnan(sig)){
            printf("  %-8.4f %-6d %-6d %-12s\n",
                   csl[i].theta_deg, csl[i].Sigma, csl[i].Sigma_true, "NaN");
            fprintf(out,"%.5f,nan,%d,%d,%d,%d,%d,%d,%d,nan,nan,%.6f,%d,%.5f\n",
                    csl[i].theta_deg, csl[i].Sigma, csl[i].Sigma_true,
                    csl[i].m, csl[i].n, csl[i].nx, csl[i].ny,
                    N, pe, conv, dmin);
            fflush(stdout); fflush(out);
            continue;
        }

        printf("  %-8.4f %-6d %-6d %-12.6f %-10.4f %-10.4f %-8d %-8.4f\n",
               csl[i].theta_deg, csl[i].Sigma, csl[i].Sigma_true,
               sig, sx, sy, N, dmin);
        if(sig < 0.0)
            printf("  WARN: negative sigma at theta=%.4f (Sig=%d)\n",
                   csl[i].theta_deg, csl[i].Sigma);

        fprintf(out,"%.5f,%.6f,%d,%d,%d,%d,%d,%d,%d,%.5f,%.5f,%.6f,%d,%.5f\n",
                csl[i].theta_deg, sig, csl[i].Sigma, csl[i].Sigma_true,
                csl[i].m, csl[i].n, csl[i].nx, csl[i].ny,
                N, sx, sy, pe, conv, dmin);
        fflush(out);
    }
    fclose(out);
    free(csl);

    printf("\n==========================================================\n");
    printf("Done.  Results -> sigma_gb_csl.csv\n");
    printf("==========================================================\n");
    return 0;
}