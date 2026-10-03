#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <omp.h>

#define PI  3.14159265358979323846

#define CUTOFF          3.0
#define OVERLAP_FRAC    0.5

#define FIRE_DT0        1e-3
#define FIRE_DTMAX      0.05
#define FIRE_DTMIN      1e-8
#define FIRE_DISP_CAP   0.05
#define FIRE_ALPHA0     0.1
#define FIRE_NMIN       5
#define FIRE_FTOL       1e-4
#define SHORT_FTOL      5e-4

#define COARSE_RELAX_STEPS  500
#define FULL_RELAX_STEPS    15000

#define SAFETY_D_PER_SIDE   8.0
#define MIN_LX_HALF_A       30.0
#define MIN_SPACINGS_Y      3.0
#define MIN_NB              2

#define A0_STAGE1       1.113895
#define ALPHA           0.0
#define T_TARGET        0.52

/* Cap on number of coarse candidates before we start coarsening the step */
#define MAX_COARSE_CANDIDATES  3000
/* Number of top coarse candidates fed into the fine sweep */
#define N_TOP_CANDIDATES       30

typedef struct { double theta_deg; int p, q; const char *label; } AngleEntry;
static const AngleEntry ANGLES[] = {
    {  4.8381, 20, 1, "very low angle   (Sigma=421)"            },
    { 10.4174,  9, 1, "low angle        (Sigma=91)"             },
    { 16.4264, 11, 2, "low-middle       (Sigma=147, S_true=49)" },
    { 21.7868,  4, 1, "upper-middle     (Sigma=21,  S_true=7)"  },
    { 25.4611, 10, 3, "plateau start    (Sigma=139)"            },
    { 32.2042,  5, 2, "deep plateau     (Sigma=39,  S_true=13)" },
};
static const int N_ANGLES = (int)(sizeof(ANGLES)/sizeof(ANGLES[0]));

typedef struct { double *x, *y; int n, cap; } PointArray;

PointArray* pa_create(int cap){
    PointArray *a = malloc(sizeof(PointArray));
    a->x = malloc(cap*sizeof(double));
    a->y = malloc(cap*sizeof(double));
    a->n = 0; a->cap = cap;
    return a;
}
void pa_free(PointArray *a){ if(a){free(a->x); free(a->y); free(a);} }
void pa_append(PointArray *a, double x, double y){
    if(a->n >= a->cap){
        a->cap = a->cap ? a->cap*2 : 8;
        a->x = realloc(a->x, a->cap*sizeof(double));
        a->y = realloc(a->y, a->cap*sizeof(double));
    }
    a->x[a->n]=x; a->y[a->n]=y; a->n++;
}
void pa_copy(PointArray *src, PointArray *dst){
    dst->n = 0;
    for(int i=0;i<src->n;i++) pa_append(dst, src->x[i], src->y[i]);
}

typedef struct {
    double ax, ay, bx, by;
    double origin_x, origin_y;
    double inv_xx, inv_xy, inv_yx, inv_yy, det;
} TriclinicBox;

void box_init(TriclinicBox *b, double ax, double ay, double bx, double by,
              double ox, double oy){
    b->ax=ax; b->ay=ay; b->bx=bx; b->by=by;
    b->origin_x=ox; b->origin_y=oy;
    double det = ax*by - ay*bx;
    b->det = det;
    b->inv_xx =  by/det;
    b->inv_xy = -bx/det;
    b->inv_yx = -ay/det;
    b->inv_yy =  ax/det;
}
void cart_to_frac(TriclinicBox *b, double x, double y, double *sx, double *sy){
    double dx = x - b->origin_x, dy = y - b->origin_y;
    *sx = b->inv_xx*dx + b->inv_xy*dy;
    *sy = b->inv_yx*dx + b->inv_yy*dy;
    *sx -= floor(*sx); *sy -= floor(*sy);
}
void frac_to_cart(TriclinicBox *b, double sx, double sy, double *x, double *y){
    *x = b->origin_x + sx*b->ax + sy*b->bx;
    *y = b->origin_y + sx*b->ay + sy*b->by;
}
void min_image(TriclinicBox *b, double xi,double yi,double xj,double yj,
               double *dx, double *dy){
    double sxi,syi,sxj,syj;
    cart_to_frac(b, xi,yi, &sxi,&syi);
    cart_to_frac(b, xj,yj, &sxj,&syj);
    double dsx = sxi-sxj, dsy = syi-syj;
    dsx -= round(dsx); dsy -= round(dsy);
    *dx = dsx*b->ax + dsy*b->bx;
    *dy = dsx*b->ay + dsy*b->by;
}

PointArray* generate_rotated_triangular_lattice_triclinic(double a, double theta,
                                                          TriclinicBox *box){
    double e1x = a,     e1y = 0.0;
    double e2x = 0.5*a, e2y = a*sqrt(3.0)/2.0;
    double c = cos(theta), s = sin(theta);
    double re1x = c*e1x - s*e1y, re1y = s*e1x + c*e1y;
    double re2x = c*e2x - s*e2y, re2y = s*e2x + c*e2y;

    double T00 = box->inv_xx*re1x + box->inv_xy*re1y;
    double T01 = box->inv_xx*re2x + box->inv_xy*re2y;
    double T10 = box->inv_yx*re1x + box->inv_yy*re1y;
    double T11 = box->inv_yx*re2x + box->inv_yy*re2y;

    double detT = T00*T11 - T01*T10;
    if (fabs(detT) < 1e-15) { fprintf(stderr,"Singular lattice transform\n"); exit(1); }
    double iT00 =  T11/detT, iT01 = -T01/detT;
    double iT10 = -T10/detT, iT11 =  T00/detT;

    double corners[4][2] = {{0,0},{1,0},{0,1},{1,1}};
    double mn_i=1e30,mx_i=-1e30,mn_j=1e30,mx_j=-1e30;
    for(int k=0;k<4;k++){
        double i = iT00*corners[k][0] + iT01*corners[k][1];
        double j = iT10*corners[k][0] + iT11*corners[k][1];
        if(i<mn_i) mn_i=i; if(i>mx_i) mx_i=i;
        if(j<mn_j) mn_j=j; if(j>mx_j) mx_j=j;
    }
    int i0=(int)floor(mn_i)-1, i1=(int)ceil(mx_i)+1;
    int j0=(int)floor(mn_j)-1, j1=(int)ceil(mx_j)+1;

    PointArray *pts = pa_create(2000);
    for(int i=i0;i<=i1;i++){
        for(int j=j0;j<=j1;j++){
            double sx = T00*i + T01*j;
            double sy = T10*i + T11*j;
            if(sx >= -1e-12 && sx < 1.0-1e-12 &&
               sy >= -1e-12 && sy < 1.0-1e-12){
                double x,y;
                frac_to_cart(box, sx, sy, &x, &y);
                pa_append(pts, x, y);
            }
        }
    }
    return pts;
}

static void check_commensurate(TriclinicBox *box, double theta_rot, double a,
                               const char *label){
    double e1x=a, e1y=0.0, e2x=0.5*a, e2y=a*sqrt(3.0)/2.0;
    double c=cos(theta_rot), s=sin(theta_rot);
    double re1x=c*e1x-s*e1y, re1y=s*e1x+c*e1y;
    double re2x=c*e2x-s*e2y, re2y=s*e2x+c*e2y;
    double det = re1x*re2y - re1y*re2x;
    double Ai = ( re2y*box->ax - re2x*box->ay)/det;
    double Aj = (-re1y*box->ax + re1x*box->ay)/det;
    double Bi = ( re2y*box->bx - re2x*box->by)/det;
    double Bj = (-re1y*box->bx + re1x*box->by)/det;
    bool A_ok = fabs(Ai-round(Ai))<1e-6 && fabs(Aj-round(Aj))<1e-6;
    bool B_ok = fabs(Bi-round(Bi))<1e-6 && fabs(Bj-round(Bj))<1e-6;
    printf("    grain %s:  A %s  B %s\n",
           label, A_ok?"yes":"NO ", B_ok?"yes":"NO ");
}

static void lj_force_energy_cell(PointArray *pos, TriclinicBox *box,
                                 double rc, double sigma, double epsilon,
                                 double *fx, double *fy, double *pe_out){
    int N = pos->n;
    double rc2 = rc*rc;
    double energy = 0.0;
    double sigma6  = pow(sigma, 6.0);
    double sigma12 = sigma6*sigma6;
    double rc6  = sigma6 / (rc2*rc2*rc2);
    double rc12 = rc6*rc6;
    double shift = 4.0*epsilon*(rc12 - rc6);
    double f_rc  = 4.0*epsilon*(12.0*rc12 - 6.0*rc6)/rc;

    memset(fx, 0, N*sizeof(double));
    memset(fy, 0, N*sizeof(double));

    double lenA = sqrt(box->ax*box->ax + box->ay*box->ay);
    double lenB = sqrt(box->bx*box->bx + box->by*box->by);
    double hA = fabs(box->det) / lenB;
    double hB = fabs(box->det) / lenA;
    int nx = (int)(hA / rc); if(nx < 1) nx = 1;
    int ny = (int)(hB / rc); if(ny < 1) ny = 1;
    int ncell = nx*ny;
    int *head = malloc(ncell*sizeof(int));
    int *lscl = malloc(N*sizeof(int));
    for(int c=0;c<ncell;c++) head[c] = -1;
    for(int i=0;i<N;i++){
        double sx,sy; cart_to_frac(box, pos->x[i], pos->y[i], &sx, &sy);
        int cx = (int)(sx*nx); int cy = (int)(sy*ny);
        if(cx<0) cx=0; if(cx>=nx) cx=nx-1;
        if(cy<0) cy=0; if(cy>=ny) cy=ny-1;
        int cell = cy*nx+cx;
        lscl[i] = head[cell]; head[cell] = i;
    }

    for(int i=0;i<N;i++){
        double sxi,syi; cart_to_frac(box, pos->x[i], pos->y[i], &sxi, &syi);
        int cx = (int)(sxi*nx); int cy = (int)(syi*ny);
        if(cx<0) cx=0; if(cx>=nx) cx=nx-1;
        if(cy<0) cy=0; if(cy>=ny) cy=ny-1;
        for(int ii=cx-1; ii<=cx+1; ii++){
            int iix = ii%nx; if(iix<0) iix += nx;
            for(int jj=cy-1; jj<=cy+1; jj++){
                int jjy = jj%ny; if(jjy<0) jjy += ny;
                int cell = jjy*nx + iix;
                for(int j=head[cell]; j!=-1; j=lscl[j]){
                    if(j <= i) continue;
                    double dx,dy;
                    min_image(box, pos->x[i],pos->y[i], pos->x[j],pos->y[j], &dx,&dy);
                    double r2 = dx*dx + dy*dy;
                    if(r2 < rc2 && r2 > 1e-12){
                        double r2i = 1.0/r2;
                        double r6  = sigma6 *r2i*r2i*r2i;
                        double r12 = r6*r6;
                        double r = sqrt(r2);
                        double f_raw = 4.0*epsilon*(12.0*r12 - 6.0*r6)/r;
                        double f_eff = f_raw - f_rc;
                        double fmag = f_eff / r;
                        double fxi = fmag*dx, fyi = fmag*dy;
                        fx[i] += fxi; fy[i] += fyi;
                        fx[j] -= fxi; fy[j] -= fyi;
                        energy += 4.0*epsilon*(r12 - r6) - shift - f_rc*(rc - r);
                    }
                }
            }
        }
    }
    free(head); free(lscl);
    *pe_out = energy;
}

static int fire_minimise(PointArray *pos, TriclinicBox *box,
                         double rc, double sigma, double epsilon,
                         int max_steps, double ftol, double *pe_out)
{
    int N = pos->n;
    double *vx = calloc(N, sizeof(double));
    double *vy = calloc(N, sizeof(double));
    double *fx = malloc(N * sizeof(double));
    double *fy = malloc(N * sizeof(double));

    double dt = FIRE_DT0, alpha = FIRE_ALPHA0;
    int n_pos = 0;
    double pe;

    lj_force_energy_cell(pos, box, rc, sigma, epsilon, fx, fy, &pe);

    int converged = 0;
    for(int step = 0; step < max_steps; step++){
        for(int i=0;i<N;i++){
            vx[i] += 0.5*dt*fx[i];
            vy[i] += 0.5*dt*fy[i];
        }
        for(int i=0;i<N;i++){
            double dx = dt*vx[i], dy = dt*vy[i];
            double d2 = dx*dx + dy*dy;
            double cap2 = FIRE_DISP_CAP*FIRE_DISP_CAP;
            if(d2 > cap2){
                double scale = sqrt(cap2/d2);
                dx *= scale; dy *= scale;
                vx[i] *= scale; vy[i] *= scale;
            }
            pos->x[i] += dx;
            pos->y[i] += dy;
        }
        for(int i=0;i<N;i++){
            double sx,sy; cart_to_frac(box, pos->x[i], pos->y[i], &sx, &sy);
            sx -= floor(sx); sy -= floor(sy);
            frac_to_cart(box, sx, sy, &pos->x[i], &pos->y[i]);
        }
        lj_force_energy_cell(pos, box, rc, sigma, epsilon, fx, fy, &pe);
        for(int i=0;i<N;i++){
            vx[i] += 0.5*dt*fx[i];
            vy[i] += 0.5*dt*fy[i];
        }
        double P=0, v2=0, f2=0;
        for(int i=0;i<N;i++){
            P  += fx[i]*vx[i] + fy[i]*vy[i];
            v2 += vx[i]*vx[i] + vy[i]*vy[i];
            f2 += fx[i]*fx[i] + fy[i]*fy[i];
        }
        double vn = sqrt(v2), fn = sqrt(f2);
        if(P > 0.0){
            n_pos++;
            if(n_pos > FIRE_NMIN){
                dt *= 1.1;
                if(dt > FIRE_DTMAX) dt = FIRE_DTMAX;
                alpha *= 0.99;
                if(alpha < 0.01) alpha = 0.01;
            }
            if(vn > 1e-14 && fn > 1e-14){
                double scale = alpha*vn/fn;
                for(int i=0;i<N;i++){
                    vx[i] = (1.0-alpha)*vx[i] + scale*fx[i];
                    vy[i] = (1.0-alpha)*vy[i] + scale*fy[i];
                }
            }
        } else {
            n_pos = 0;
            dt *= 0.5;
            if(dt < FIRE_DTMIN) dt = FIRE_DTMIN;
            alpha = FIRE_ALPHA0;
            memset(vx, 0, N*sizeof(double));
            memset(vy, 0, N*sizeof(double));
        }
        if((step & 127) == 0){
            double fmax2 = 0;
            for(int i=0;i<N;i++){
                double ff = fx[i]*fx[i] + fy[i]*fy[i];
                if(ff > fmax2) fmax2 = ff;
            }
            if(sqrt(fmax2) < ftol){ converged = 1; break; }
        }
    }
    if(!converged){
        double fmax2 = 0;
        for(int i=0;i<N;i++){
            double ff = fx[i]*fx[i] + fy[i]*fy[i];
            if(ff > fmax2) fmax2 = ff;
        }
        if(sqrt(fmax2) < ftol) converged = 1;
    }
    if(pe_out) *pe_out = pe;
    free(vx); free(vy); free(fx); free(fy);
    return converged;
}

void remove_overlaps(PointArray *A, double rcut, TriclinicBox *box){
    int N = A->n;
    if(N < 2) return;
    double rcut2 = rcut*rcut;
    bool *keep = malloc(N*sizeof(bool));
    for(int i=0;i<N;i++) keep[i] = true;

    for(int i=0;i<N;i++){
        if(!keep[i]) continue;
        for(int j=i+1;j<N;j++){
            if(!keep[j]) continue;
            double dx,dy;
            min_image(box, A->x[i],A->y[i], A->x[j],A->y[j], &dx,&dy);
            if(dx*dx + dy*dy < rcut2){
                if( ((i + j) & 1) == 0 ) keep[j] = false;
                else                       keep[i] = false;
            }
        }
    }
    int Nn = 0;
    for(int i=0;i<N;i++) if(keep[i]){
        A->x[Nn] = A->x[i];
        A->y[Nn] = A->y[i];
        Nn++;
    }
    A->n = Nn;
    free(keep);
}

double compute_bulk_energy_per_atom(double a, double theta, TriclinicBox *box,
                                    double rc, double sigma, double epsilon){
    PointArray *bulk = generate_rotated_triangular_lattice_triclinic(a, theta, box);
    if(bulk->n == 0){ pa_free(bulk); return 0.0; }
    remove_overlaps(bulk, 0.01*a, box);
    double pe;
    fire_minimise(bulk, box, rc, sigma, epsilon, 3000, 1e-4, &pe);
    double e = pe / bulk->n;
    pa_free(bulk);
    return e;
}

typedef struct { PointArray *relaxed; double energy; } TranslationResult;

TranslationResult* process_translation(double dx, double dy,
                                       PointArray *left, PointArray *right,
                                       TriclinicBox *box,
                                       double rc, double sigma, double epsilon, double a){
    TranslationResult *res = malloc(sizeof(TranslationResult));
    double hdx = 0.5*dx, hdy = 0.5*dy;

    PointArray *atoms = pa_create(left->n + right->n);
    for(int i=0;i<left->n;i++){
        double x = left->x[i] - hdx, y = left->y[i] - hdy;
        double sx,sy; cart_to_frac(box, x,y, &sx,&sy);
        sx -= floor(sx); sy -= floor(sy);
        frac_to_cart(box, sx, sy, &x, &y);
        pa_append(atoms, x, y);
    }
    for(int i=0;i<right->n;i++){
        double x = right->x[i] + hdx, y = right->y[i] + hdy;
        double sx,sy; cart_to_frac(box, x,y, &sx,&sy);
        sx -= floor(sx); sy -= floor(sy);
        frac_to_cart(box, sx, sy, &x, &y);
        pa_append(atoms, x, y);
    }
    if(atoms->n == 0){
        res->relaxed = NULL; res->energy = 1e30;
        pa_free(atoms);
        return res;
    }
    remove_overlaps(atoms, OVERLAP_FRAC*a, box);

    double pe;
    fire_minimise(atoms, box, rc, sigma, epsilon,
                  COARSE_RELAX_STEPS, SHORT_FTOL, &pe);

    res->relaxed = atoms;
    res->energy  = pe;
    return res;
}

/* ============================================================ */
int main(int argc, char **argv){
    int angle_idx = 1;
    if(argc > 1) angle_idx = atoi(argv[1]);
    if(angle_idx < 1 || angle_idx > N_ANGLES){
        fprintf(stderr, "Usage: %s [1..%d]\n", argv[0], N_ANGLES);
        for(int i=0;i<N_ANGLES;i++)
            fprintf(stderr,"  %d : %6.4f deg  %s\n", i+1,
                    ANGLES[i].theta_deg, ANGLES[i].label);
        return 1;
    }

    double sigma = 1.0, epsilon = 1.0;
    double a0 = A0_STAGE1;
    double a  = a0 * (1.0 + ALPHA*T_TARGET);

    const AngleEntry *AE = &ANGLES[angle_idx-1];
    int p = AE->p, q = AE->q;
    double theta_deg = AE->theta_deg;
    double theta = 2.0*atan2(sqrt(3.0)*q, 2.0*p + q);
    double half_theta = 0.5*theta;
    int sigma_val  = p*p + p*q + q*q;
    int sigma_true = (sigma_val % 3 == 0) ? sigma_val/3 : sigma_val;

    printf("=== Bicrystal run %d/%d ===\n", angle_idx, N_ANGLES);
    printf("  label      : %s\n", AE->label);
    printf("  theta      : %.4f deg  (p=%d, q=%d)\n", theta_deg, p, q);
    printf("  Sigma      : %d   (Sigma_true = %d)\n", sigma_val, sigma_true);
    printf("  a          : %.6f\n", a);

    double theta_rad = theta;
    double D_over_a  = 1.0 / theta_rad;
    double Lx_half_target_a = SAFETY_D_PER_SIDE * D_over_a;
    if(Lx_half_target_a < MIN_LX_HALF_A) Lx_half_target_a = MIN_LX_HALF_A;
    double Lx_target_a = 2.0 * Lx_half_target_a;

    double Lx_base_a = sqrt((double)sigma_val);
    double Ly_base_a = sqrt(3.0*(double)sigma_val);

    int N_A = (int)ceil(Lx_target_a / Lx_base_a);
    if(N_A < 2) N_A = 2;
    if(N_A % 2 != 0) N_A++;

    double spacings_per_Ly_base = Ly_base_a / D_over_a;
    int N_B = (int)ceil(MIN_SPACINGS_Y / spacings_per_Ly_base);
    if(N_B < MIN_NB) N_B = MIN_NB;
    if(N_B % 2 != 0) N_B++;

    double box_ax = N_A * Lx_base_a * a;
    double box_ay = 0.0;
    double box_bx = 0.0;
    double box_by = N_B * Ly_base_a * a;
    double origin_x = -0.5*(box_ax + box_bx);
    double origin_y = -0.5*(box_ay + box_by);

    TriclinicBox box;
    box_init(&box, box_ax, box_ay, box_bx, box_by, origin_x, origin_y);

    double Lx_half_a = 0.5*box_ax / a;
    double Ly_a      = box_by / a;
    int    N_atoms_expected = 2 * N_A * N_B * sigma_val;

    printf("\n  Adaptive box:\n");
    printf("    dislocation spacing D = %.3f a\n", D_over_a);
    printf("    N_A = %d   (Lx = %d * sqrt(Sigma) * a)\n", N_A, N_A);
    printf("    N_B = %d   (Ly = %d * sqrt(3 Sigma) * a)\n", N_B, N_B);
    printf("    Lx = %.2f a   Ly = %.2f a\n", box_ax/a, Ly_a);
    printf("    Lx/2 = %.2f a  ->  Lx/(2D) = %.2f\n",
           Lx_half_a, Lx_half_a / D_over_a);
    printf("    expected N atoms = %d\n", N_atoms_expected);
    printf("  Commensurability:\n");
    check_commensurate(&box, +half_theta, a, "1 (+theta/2)");
    check_commensurate(&box, -half_theta, a, "2 (-theta/2)");

    double rc = CUTOFF*sigma;

    printf("\nComputing bulk reference (FIRE)...\n");
    double e_plus  = compute_bulk_energy_per_atom(a, +half_theta, &box, rc, sigma, epsilon);
    double e_minus = compute_bulk_energy_per_atom(a, -half_theta, &box, rc, sigma, epsilon);
    double e_bulk  = 0.5*(e_plus + e_minus);
    printf("  e_bulk = %.6f  (Stage-1 e0 = -3.147364)\n", e_bulk);

    PointArray *all_L = generate_rotated_triangular_lattice_triclinic(a, +half_theta, &box);
    PointArray *all_R = generate_rotated_triangular_lattice_triclinic(a, -half_theta, &box);
    PointArray *left  = pa_create(all_L->n);
    PointArray *right = pa_create(all_R->n);
    for(int i=0;i<all_L->n;i++){
        double sx,sy; cart_to_frac(&box, all_L->x[i], all_L->y[i], &sx, &sy);
        if(sx < 0.5) pa_append(left, all_L->x[i], all_L->y[i]);
    }
    for(int i=0;i<all_R->n;i++){
        double sx,sy; cart_to_frac(&box, all_R->x[i], all_R->y[i], &sx, &sy);
        if(sx >= 0.5) pa_append(right, all_R->x[i], all_R->y[i]);
    }
    pa_free(all_L); pa_free(all_R);
    printf("\nAtoms before overlap removal: left=%d right=%d total=%d\n",
           left->n, right->n, left->n + right->n);

    double interface_length = box_by;

    /* ============================================================
       COARSE SWEEP — covers one full CSL period in x and y.
       The CSL periods along the box axes are:
          P_x = sqrt(Sigma)   * a
          P_y = sqrt(3 Sigma) * a
       We sweep  dx in [-P_x/2, +P_x/2] and dy in [-P_y/2, +P_y/2],
       centered on zero.  If the resulting grid exceeds
       MAX_COARSE_CANDIDATES, the step is increased adaptively.
       ============================================================ */
    double csl_period_x = sqrt((double)sigma_val) * a;
    double csl_period_y = sqrt(3.0*(double)sigma_val) * a;

    double coarse_step = 0.4 * a;
    int coarse_n_dx = (int)ceil(0.5 * csl_period_x / coarse_step);
    int coarse_n_dy = (int)ceil(0.5 * csl_period_y / coarse_step);
    int coarse_max  = (2*coarse_n_dx + 1) * (2*coarse_n_dy + 1);

    while(coarse_max > MAX_COARSE_CANDIDATES){
        coarse_step *= 1.15;
        coarse_n_dx = (int)ceil(0.5 * csl_period_x / coarse_step);
        coarse_n_dy = (int)ceil(0.5 * csl_period_y / coarse_step);
        coarse_max  = (2*coarse_n_dx + 1) * (2*coarse_n_dy + 1);
    }

    printf("\nCoarse sweep (full CSL period):\n");
    printf("  CSL period x = %.3f a, y = %.3f a\n",
           csl_period_x/a, csl_period_y/a);
    printf("  step = %.3f a  ->  grid %d x %d = %d candidates\n",
           coarse_step/a, 2*coarse_n_dx+1, 2*coarse_n_dy+1, coarse_max);

    typedef struct { double dx, dy, gb; } Cand;
    Cand *cands = malloc(coarse_max * sizeof(Cand));
    int ncand = 0;

    #pragma omp parallel for collapse(2) schedule(dynamic)
    for(int ix = -coarse_n_dx; ix <= coarse_n_dx; ix++){
        for(int iy = -coarse_n_dy; iy <= coarse_n_dy; iy++){
            double dx = ix*coarse_step;
            double dy = iy*coarse_step;
            TranslationResult *r = process_translation(dx, dy, left, right,
                                                       &box, rc, sigma, epsilon, a);
            if(!r || !r->relaxed){ if(r) free(r); continue; }
            double gb = (r->energy - r->relaxed->n * e_bulk) / (2.0 * interface_length);
            #pragma omp critical
            {
                if(gb >= 0.0 && ncand < coarse_max){
                    cands[ncand].dx = dx;
                    cands[ncand].dy = dy;
                    cands[ncand].gb = gb;
                    ncand++;
                }
            }
            pa_free(r->relaxed); free(r);
        }
    }
    printf("  kept %d non-negative candidates\n", ncand);
    if(ncand == 0){
        printf("No valid candidates. Exiting.\n");
        free(cands); pa_free(left); pa_free(right);
        return 1;
    }
    for(int i=0;i<ncand-1;i++)
        for(int j=i+1;j<ncand;j++)
            if(cands[j].gb < cands[i].gb){
                Cand t = cands[i]; cands[i] = cands[j]; cands[j] = t;
            }

    int top = N_TOP_CANDIDATES;
    if(top > ncand) top = ncand;
    printf("Top coarse candidates:\n");
    for(int i=0;i<top && i<ncand;i++)
        printf("  %d  dx=%8.4f  dy=%8.4f  GB=%10.6f\n",
               i+1, cands[i].dx/a, cands[i].dy/a, cands[i].gb);

    /* ============================================================
       FINE SWEEP — step 0.08a, half-range 0.5a around each of the
       top N_TOP_CANDIDATES coarse points.  Guaranteed to cover the
       gap left by the coarse step.
       ============================================================ */
    double fine_step  = 0.08*a;
    double fine_delta = 0.5*a;
    int fine_n = (int)(fine_delta / fine_step);

    double best_gb = 1e30, best_dx = 0, best_dy = 0;
    PointArray *best_relaxed = NULL;

    for(int c=0; c<top && c<ncand; c++){
        double cdx = cands[c].dx, cdy = cands[c].dy;
        printf("\nFine sweep around candidate %d/%d (dx=%.4f, dy=%.4f)\n",
               c+1, top, cdx/a, cdy/a);
        #pragma omp parallel for collapse(2) schedule(dynamic)
        for(int ix = -fine_n; ix <= fine_n; ix++){
            for(int iy = -fine_n; iy <= fine_n; iy++){
                double dx = cdx + ix*fine_step;
                double dy = cdy + iy*fine_step;
                TranslationResult *r = process_translation(dx, dy, left, right,
                                                           &box, rc, sigma, epsilon, a);
                if(!r || !r->relaxed){ if(r) free(r); continue; }
                double gb = (r->energy - r->relaxed->n * e_bulk) / (2.0 * interface_length);
                #pragma omp critical
                {
                    if(gb < best_gb && gb >= 0.0){
                        best_gb = gb; best_dx = dx; best_dy = dy;
                        if(best_relaxed) pa_free(best_relaxed);
                        best_relaxed = pa_create(r->relaxed->n);
                        pa_copy(r->relaxed, best_relaxed);
                        printf("  new best: GB=%.6f  dx=%.4f  dy=%.4f\n",
                               gb, dx/a, dy/a);
                    }
                }
                pa_free(r->relaxed); free(r);
            }
        }
    }
    free(cands);

    if(!best_relaxed){
        fprintf(stderr, "ERROR: no valid fine-sweep candidate found. Exiting.\n");
        pa_free(left); pa_free(right);
        return 1;
    }

    printf("\nFull FIRE relaxation on best candidate...\n");
    double pe_final;
    fire_minimise(best_relaxed, &box, rc, sigma, epsilon,
                  FULL_RELAX_STEPS, FIRE_FTOL, &pe_final);

    double sigma_gb = (pe_final - best_relaxed->n * e_bulk) / (2.0 * interface_length);

    printf("\n=== Final result ===\n");
    printf("  theta = %.4f deg  (Sigma=%d, Sigma_true=%d)\n",
           theta_deg, sigma_val, sigma_true);
    printf("  N atoms      : %d\n", best_relaxed->n);
    printf("  Lx / Ly      : %.2f a / %.2f a\n", box_ax/a, box_by/a);
    printf("  Lx/(2D)      : %.2f  (SAFETY, want > %.1f)\n",
           0.5*box_ax/(a*D_over_a), SAFETY_D_PER_SIDE);
    printf("  best shift   : dx = %.4f a, dy = %.4f a\n",
           best_dx/a, best_dy/a);
    printf("  sigma_GB     : %.6f\n", sigma_gb);

    char fname[256];
    snprintf(fname, sizeof(fname),
             "bicrystal_Sigma%d_theta%.4fdeg.dat", sigma_val, theta_deg);
    FILE *fp = fopen(fname, "w");
    if(fp){
        fprintf(fp,"# Bicrystal: theta = %.4f deg, Sigma = %d, Sigma_true = %d\n",
                theta_deg, sigma_val, sigma_true);
        fprintf(fp,"# label    = %s\n", AE->label);
        fprintf(fp,"# a        = %.12f\n", a);
        fprintf(fp,"# box_ax   = %.12f  box_ay = %.12f\n", box.ax, box.ay);
        fprintf(fp,"# box_bx   = %.12f  box_by = %.12f\n", box.bx, box.by);
        fprintf(fp,"# origin_x = %.12f  origin_y = %.12f\n",
                box.origin_x, box.origin_y);
        fprintf(fp,"# N_atoms  = %d\n", best_relaxed->n);
        fprintf(fp,"# sigma_GB = %.12f (per interface)\n", sigma_gb);
        fprintf(fp,"# dx = %.6f a   dy = %.6f a\n", best_dx/a, best_dy/a);
        fprintf(fp,"# Columns: x y grain\n");
        for(int i=0;i<best_relaxed->n;i++){
            double sx,sy;
            cart_to_frac(&box, best_relaxed->x[i], best_relaxed->y[i], &sx, &sy);
            int grain = (sx < 0.5) ? 1 : 2;
            fprintf(fp,"%.12f %.12f %d\n",
                    best_relaxed->x[i], best_relaxed->y[i], grain);
        }
        fclose(fp);
        printf("\n  Structure written to %s\n", fname);
    }
    pa_free(best_relaxed);
    pa_free(left);
    pa_free(right);
    return 0;
}