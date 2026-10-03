/* ============================================================
   melt_scan.c -- premelting study version (v6.5, OMP over T)

   Changes vs v6.4:
     - identify_bins_from_profile() now smooths the per-bin
       psi6 and rho profiles with a 5-bin moving average before
       searching for the min (GB) and max (bulk) bins. This
       removes the noise-driven bin selection that caused the
       "GB more ordered than bulk" anomaly.
     - No other changes. The MSD block remains removed.

   Compile:
     gcc -O2 -fopenmp -o melt_scan melt_scan.c -lm

   Run:
     OMP_NUM_THREADS=4 ./melt_scan perfect_Sigma39_theta32.2042deg.dat \
                                   bicrystal_Sigma39_theta32.2042deg.dat \
                                   0.28 0.44 0.02 500000 1.0 0.005 5 0.0 1
   ============================================================ */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <omp.h>

#define PI 3.14159265358979323846

#define CUTOFF      3.0
#define KB          1.0

#define DT_DEFAULT       0.005
#define GAMMA_DEFAULT    1.0
#define SAMPLE_EVERY     100
#define SAMPLE_COORD     500
#define N_BINS           60

#define BARO_INTERVAL       20
#define BARO_KAPPA_DEFAULT  1e-2
#define BARO_MAX_SCALE      0.02

#define NN_CUTOFF_FRAC      1.35

#define GB_WIN_HALF     2
#define BULK_WIN_HALF   2

#define MIN_BIN_SEP     10
#define RHO_MIN_SMOOTH  0.3

/* Smoothing half-width for profile analysis (5-bin moving average) */
#define SMOOTH_HALF     2

/* ---------- thread-local RNG ---------- */
static __thread unsigned long long g_rng_state = 0x853c49e6748fea9bULL;
static void rng_seed(unsigned long long s){ g_rng_state = s ? s : 0x853c49e6748fea9bULL; }
static inline double rng_uniform(void){
    unsigned long long x = g_rng_state;
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    g_rng_state = x;
    return (double)((x * 2685821657736338717ULL) >> 11) * (1.0/9007199254740992.0);
}
static inline double rng_normal(void){
    double u1 = rng_uniform(); if(u1 < 1e-300) u1 = 1e-300;
    double u2 = rng_uniform();
    return sqrt(-2.0*log(u1)) * cos(2.0*PI*u2);
}

/* ---------- PointArray ---------- */
typedef struct { double *x, *y; int n, cap; } PointArray;
static PointArray* pa_create(int cap){
    PointArray *a = malloc(sizeof(PointArray));
    if(cap < 1) cap = 1;
    a->x = malloc(cap * sizeof(double));
    a->y = malloc(cap * sizeof(double));
    a->n = 0; a->cap = cap;
    return a;
}
static void pa_free(PointArray *a){ if(a){free(a->x); free(a->y); free(a);} }
static void pa_append(PointArray *a, double x, double y){
    if(a->n >= a->cap){
        a->cap *= 2;
        a->x = realloc(a->x, a->cap*sizeof(double));
        a->y = realloc(a->y, a->cap*sizeof(double));
    }
    a->x[a->n]=x; a->y[a->n]=y; a->n++;
}
static void pa_copy(const PointArray *src, PointArray *dst){
    dst->n = 0;
    for(int i=0;i<src->n;i++) pa_append(dst, src->x[i], src->y[i]);
}

/* ---------- triclinic box ---------- */
typedef struct {
    double ax, ay, bx, by;
    double origin_x, origin_y;
    double inv_xx, inv_xy, inv_yx, inv_yy, det;
    double Lx, Ly;
} TriclinicBox;

static void box_init(TriclinicBox *b, double ax, double ay, double bx, double by,
                     double ox, double oy){
    b->ax=ax; b->ay=ay; b->bx=bx; b->by=by;
    b->origin_x=ox; b->origin_y=oy;
    double det = ax*by - ay*bx;
    b->det = det;
    b->inv_xx =  by/det;
    b->inv_xy = -bx/det;
    b->inv_yx = -ay/det;
    b->inv_yy =  ax/det;
    b->Lx = sqrt(ax*ax + ay*ay);
    b->Ly = sqrt(bx*bx + by*by);
}
static void cart_to_frac(const TriclinicBox *b, double x, double y, double *sx, double *sy){
    double dx = x - b->origin_x, dy = y - b->origin_y;
    *sx = b->inv_xx*dx + b->inv_xy*dy;
    *sy = b->inv_yx*dx + b->inv_yy*dy;
    *sx -= floor(*sx); *sy -= floor(*sy);
}
static void frac_to_cart(const TriclinicBox *b, double sx, double sy, double *x, double *y){
    *x = b->origin_x + sx*b->ax + sy*b->bx;
    *y = b->origin_y + sx*b->ay + sy*b->by;
}
static void min_image(const TriclinicBox *b, double xi,double yi,double xj,double yj,
                      double *dx, double *dy){
    double sxi,syi,sxj,syj;
    cart_to_frac(b, xi,yi, &sxi,&syi);
    cart_to_frac(b, xj,yj, &sxj,&syj);
    double dsx = sxi-sxj, dsy = syi-syj;
    dsx -= round(dsx); dsy -= round(dsy);
    *dx = dsx*b->ax + dsy*b->bx;
    *dy = dsx*b->ay + dsy*b->by;
}
static void wrap_pbc(const TriclinicBox *b, double *x, double *y){
    double sx, sy;
    cart_to_frac(b, *x, *y, &sx, &sy);
    frac_to_cart(b, sx, sy, x, y);
}
static void rescale_box(TriclinicBox *box, PointArray *pos, double s){
    double cx = box->origin_x + 0.5*(box->ax + box->bx);
    double cy = box->origin_y + 0.5*(box->ay + box->by);
    double new_ax = box->ax * s, new_ay = box->ay * s;
    double new_bx = box->bx * s, new_by = box->by * s;
    double new_ox = cx - 0.5*(new_ax + new_bx);
    double new_oy = cy - 0.5*(new_ay + new_by);
    for(int i=0;i<pos->n;i++){
        pos->x[i] = cx + s * (pos->x[i] - cx);
        pos->y[i] = cy + s * (pos->y[i] - cy);
    }
    box_init(box, new_ax, new_ay, new_bx, new_by, new_ox, new_oy);
}

/* ---------- force + energy + virial ---------- */
static double compute_force_energy(PointArray *pos, const TriclinicBox *box,
                                   double rc, double sigma, double epsilon,
                                   double *fx, double *fy,
                                   int *head, int *lscl, int nx, int ny,
                                   double *virial_out)
{
    int N = pos->n;
    double rc2 = rc*rc;
    double sigma6  = pow(sigma, 6.0);
    double sigma12 = sigma6*sigma6;
    double rc6  = sigma6 / (rc2*rc2*rc2);
    double rc12 = rc6*rc6;
    double shift = 4.0*epsilon*(rc12 - rc6);
    double f_rc  = 4.0*epsilon*(12.0*rc12 - 6.0*rc6)/rc;

    memset(fx, 0, N*sizeof(double));
    memset(fy, 0, N*sizeof(double));

    int ncell = nx*ny;
    for(int c=0;c<ncell;c++) head[c] = -1;
    for(int i=0;i<N;i++){
        double sx, sy;
        cart_to_frac(box, pos->x[i], pos->y[i], &sx, &sy);
        int cx = (int)(sx*nx); int cy = (int)(sy*ny);
        if(cx<0) cx=0; if(cx>=nx) cx=nx-1;
        if(cy<0) cy=0; if(cy>=ny) cy=ny-1;
        int cell = cy*nx + cx;
        lscl[i] = head[cell]; head[cell] = i;
    }

    double energy = 0.0;
    double virial = 0.0;

    for(int i=0;i<N;i++){
        double sxi, syi;
        cart_to_frac(box, pos->x[i], pos->y[i], &sxi, &syi);
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
                    double dx, dy;
                    min_image(box, pos->x[i],pos->y[i], pos->x[j],pos->y[j], &dx,&dy);
                    double r2 = dx*dx + dy*dy;
                    if(r2 < rc2 && r2 > 1e-12){
                        double r2i = 1.0/r2;
                        double r6  = sigma6 * r2i*r2i*r2i;
                        double r12 = r6*r6;
                        double r = sqrt(r2);
                        double f_raw = 4.0*epsilon*(12.0*r12 - 6.0*r6)/r;
                        double f_eff = f_raw - f_rc;
                        double fmag = f_eff / r;
                        double fxi = fmag*dx, fyi = fmag*dy;
                        fx[i] += fxi; fy[i] += fyi;
                        fx[j] -= fxi; fy[j] -= fyi;
                        energy += 4.0*epsilon*(r12 - r6) - shift - f_rc*(rc - r);
                        virial += f_eff * r;
                    }
                }
            }
        }
    }
    if(virial_out) *virial_out = virial;
    return energy;
}

static double compute_pressure(const TriclinicBox *box,
                               const double *vx, const double *vy, int N,
                               double virial)
{
    double area = fabs(box->det);
    double KE = 0.0;
    for(int i=0;i<N;i++) KE += 0.5*(vx[i]*vx[i] + vy[i]*vy[i]);
    return (KE + 0.5 * virial) / area;
}

/* ---------- layer-resolved coord / psi4 / psi6 ---------- */
static void accumulate_layer_profile(PointArray *pos, const TriclinicBox *box,
                                     double r_nn,
                                     int *head, int *lscl, int nx, int ny,
                                     double *rho_sum,
                                     double *psi6_sum,
                                     double *psi4_sum,
                                     double *coord_sum)
{
    int N = pos->n;
    double r2cut = r_nn*r_nn;
    int ncell = nx*ny;
    for(int c=0;c<ncell;c++) head[c] = -1;
    for(int i=0;i<N;i++){
        double sx, sy;
        cart_to_frac(box, pos->x[i], pos->y[i], &sx, &sy);
        int cx = (int)(sx*nx); int cy = (int)(sy*ny);
        if(cx<0) cx=0; if(cx>=nx) cx=nx-1;
        if(cy<0) cy=0; if(cy>=ny) cy=ny-1;
        int cell = cy*nx + cx;
        lscl[i] = head[cell]; head[cell] = i;
    }

    for(int i=0;i<N;i++){
        double sxi, syi;
        cart_to_frac(box, pos->x[i], pos->y[i], &sxi, &syi);
        int b = (int)(sxi * N_BINS);
        if(b < 0) b = 0;
        if(b >= N_BINS) b = N_BINS - 1;

        rho_sum[b] += 1.0;

        int cx = (int)(sxi*nx); int cy = (int)(syi*ny);
        if(cx<0) cx=0; if(cx>=nx) cx=nx-1;
        if(cy<0) cy=0; if(cy>=ny) cy=ny-1;

        int n_nb = 0;
        double s6c = 0.0, s6s = 0.0, s4c = 0.0, s4s = 0.0;

        for(int ii=cx-1; ii<=cx+1; ii++){
            int iix = ii%nx; if(iix<0) iix += nx;
            for(int jj=cy-1; jj<=cy+1; jj++){
                int jjy = jj%ny; if(jjy<0) jjy += ny;
                int cell = jjy*nx + iix;
                for(int j=head[cell]; j!=-1; j=lscl[j]){
                    if(j == i) continue;
                    double dx, dy;
                    min_image(box, pos->x[i],pos->y[i], pos->x[j],pos->y[j], &dx,&dy);
                    double r2 = dx*dx + dy*dy;
                    if(r2 < r2cut){
                        double theta = atan2(dy, dx);
                        s6c += cos(6.0*theta); s6s += sin(6.0*theta);
                        s4c += cos(4.0*theta); s4s += sin(4.0*theta);
                        n_nb++;
                    }
                }
            }
        }
        if(n_nb > 0){
            psi6_sum[b] += sqrt(s6c*s6c + s6s*s6s) / n_nb;
            psi4_sum[b] += sqrt(s4c*s4c + s4s*s4s) / n_nb;
        }
        coord_sum[b] += n_nb;
    }
}

/* ---------- profile smoothing (circular moving average) ---------- */
static void smooth_profile(const double *in, double *out, int n, int half)
{
    for(int b = 0; b < n; b++){
        double sum = 0.0;
        int cnt = 0;
        for(int d = -half; d <= half; d++){
            int bb = (b + d + n) % n;
            sum += in[bb];
            cnt++;
        }
        out[b] = sum / cnt;
    }
}

/* ---------- Diagnostic: identify GB and bulk bins from initial structure ---------- */
static void identify_gb_and_bulk_bins(const PointArray *init,
                                      const TriclinicBox *box,
                                      int nx, int ny,
                                      int *gb_bin_out, int *bulk_bin_out)
{
    int N = init->n;
    int ncell = nx*ny;
    int    *head  = malloc(ncell * sizeof(int));
    int    *lscl  = malloc(N * sizeof(int));
    double *rho   = calloc(N_BINS, sizeof(double));
    double *psi6  = calloc(N_BINS, sizeof(double));
    double *psi4  = calloc(N_BINS, sizeof(double));
    double *coord = calloc(N_BINS, sizeof(double));

    PointArray *pos = pa_create(N);
    pa_copy(init, pos);

    double a_nn = sqrt(2.0 * box->Lx * box->Ly / (sqrt(3.0) * N));
    double r_nn = NN_CUTOFF_FRAC * a_nn;

    accumulate_layer_profile(pos, box, r_nn, head, lscl, nx, ny,
                             rho, psi6, psi4, coord);

    /* smooth before searching */
    double psi6_s[N_BINS];
    double rho_s [N_BINS];
    smooth_profile(psi6, psi6_s, N_BINS, SMOOTH_HALF);
    smooth_profile(rho,  rho_s,  N_BINS, SMOOTH_HALF);

    int gb = -1;
    double p6min = 1e30;
    for(int b = 2; b < N_BINS - 2; b++){
        if(rho_s[b] < 1.0) continue;
        if(psi6_s[b] < p6min){ p6min = psi6_s[b]; gb = b; }
    }
    if(gb < 0) gb = N_BINS / 2;

    int bulk = -1;
    double p6max = -1.0;
    for(int b = 2; b < N_BINS - 2; b++){
        if(rho_s[b] < 1.0) continue;
        int d = abs(b - gb);
        int dwrap = (d < N_BINS - d) ? d : N_BINS - d;
        if(dwrap < MIN_BIN_SEP) continue;
        if(psi6_s[b] > p6max){ p6max = psi6_s[b]; bulk = b; }
    }
    if(bulk < 0) bulk = (gb + N_BINS / 2) % N_BINS;

    *gb_bin_out   = gb;
    *bulk_bin_out = bulk;

    printf("[diag] initial-structure bin id (smoothed):\n");
    printf("[diag]   gb_bin   = %2d  (psi6 = %.4f)\n", gb, p6min);
    printf("[diag]   bulk_bin = %2d  (psi6 = %.4f)\n", bulk, p6max);
    fflush(stdout);

    free(head); free(lscl);
    free(rho); free(psi6); free(psi4); free(coord);
    pa_free(pos);
}

/* ---------- windowed average over a bin range ---------- */
static double windowed_average(const double *values, const double *counts,
                               int n_bins, int center_bin, int half_width)
{
    double sum = 0.0;
    int    cnt = 0;
    for(int db = -half_width; db <= half_width; db++){
        int b = (center_bin + db + n_bins) % n_bins;
        if(counts[b] > 0.0){
            sum += values[b];
            cnt++;
        }
    }
    return cnt > 0 ? sum / cnt : 0.0;
}

/* ---------- header parsing ---------- */
static const char *find_key(const char *line, const char *key){
    size_t klen = strlen(key);
    if(klen == 0) return NULL;
    const char *p = line;
    while((p = strstr(p, key)) != NULL){
        unsigned char prev = (p == line) ? ' ' : (unsigned char)p[-1];
        unsigned char next = (unsigned char)p[klen];
        int prev_ok = !( (prev>='a'&&prev<='z') ||
                         (prev>='A'&&prev<='Z') ||
                         (prev>='0'&&prev<='9') || prev=='_' );
        int next_ok = !( (next>='a'&&next<='z') ||
                         (next>='A'&&next<='Z') ||
                         (next>='0'&&next<='9') || next=='_' );
        if(prev_ok && next_ok) return p;
        p++;
    }
    return NULL;
}
static int get_double_after(const char *line, const char *key, double *out){
    const char *p = find_key(line, key);
    if(!p) return 0;
    p = strchr(p, '=');
    if(!p) return 0;
    char *end;
    double v = strtod(p+1, &end);
    if(end == p+1) return 0;
    *out = v;
    return 1;
}
static int get_int_after(const char *line, const char *key, int *out){
    const char *p = find_key(line, key);
    if(!p) return 0;
    p = strchr(p, '=');
    if(!p) return 0;
    char *end;
    long v = strtol(p+1, &end, 10);
    if(end == p+1) return 0;
    *out = (int)v;
    return 1;
}

static int read_dat(const char *fname, PointArray *pos,
                    double *box_ax, double *box_ay,
                    double *box_bx, double *box_by,
                    double *a_out, int *N_atoms)
{
    FILE *fp = fopen(fname, "r");
    if(!fp){ fprintf(stderr,"Cannot open %s\n", fname); return 1; }

    char line[512];
    *box_ax = *box_ay = *box_bx = *box_by = 0.0;
    *a_out = 0.0;
    *N_atoms = 0;
    int got_ax=0, got_bx=0, got_a=0, got_N=0;

    while(fgets(line, sizeof(line), fp)){
        if(line[0] != '#') break;
        double v;
        if(!got_ax && get_double_after(line, "box_ax", &v)){ *box_ax = v; got_ax = 1; }
        if(get_double_after(line, "box_ay", &v)) *box_ay = v;
        if(!got_bx && get_double_after(line, "box_bx", &v)){ *box_bx = v; got_bx = 1; }
        if(get_double_after(line, "box_by", &v)) *box_by = v;
        if(!got_a  && get_double_after(line, "a", &v)){ *a_out = v; got_a = 1; }
        int iv;
        if(!got_N && get_int_after(line, "N_atoms", &iv)){ *N_atoms = iv; got_N = 1; }
    }
    if(*box_ax == 0.0 || *box_by == 0.0){
        fprintf(stderr,"Missing box vectors in %s\n", fname);
        fclose(fp); return 1;
    }
    if(*a_out == 0.0){
        fprintf(stderr,"Missing a in %s\n", fname);
        fclose(fp); return 1;
    }
    rewind(fp);
    while(fgets(line, sizeof(line), fp)){
        if(line[0] == '#') continue;
        double x, y;
        int g;
        if(sscanf(line, "%lf %lf %d", &x, &y, &g) == 3) pa_append(pos, x, y);
        else if(sscanf(line, "%lf %lf", &x, &y) == 2)   pa_append(pos, x, y);
    }
    fclose(fp);
    if(*N_atoms == 0) *N_atoms = pos->n;
    return 0;
}

static void basename_no_ext(const char *path, char *out, size_t out_sz){
    const char *base = strrchr(path, '/');
    if(!base) base = strrchr(path, '\\');
    base = base ? base + 1 : path;
    size_t n = strlen(base);
    if(n >= out_sz) n = out_sz - 1;
    memcpy(out, base, n); out[n] = '\0';
    char *dot = strrchr(out, '.');
    if(dot) *dot = '\0';
}

typedef struct {
    double coord_avg, psi6_avg, psi4_avg;
    double pe_per_atom, e_per_atom;
    double P_avg;
    double a_avg;
    double psi6_bulk, psi6_gb;
    double psi4_bulk, psi4_gb;
    double coord_bulk, coord_gb;
    double rho_bulk, rho_gb;
    double layer_rho[N_BINS];
    double layer_psi6[N_BINS];
    double layer_psi4[N_BINS];
    double layer_coord[N_BINS];
    int n_layer_samples;
    int gb_bin, bulk_bin;
} MDResult;

/* ============================================================
   Post-processing: identify GB and bulk bins from the
   production-averaged profile of the bicrystal, using a
   smoothed version of the profile to avoid noise-driven
   bin selection.
   ============================================================ */
static void identify_bins_from_profile(const MDResult *r,
                                       int *gb_out, int *bulk_out)
{
    double psi6_s[N_BINS];
    double rho_s [N_BINS];
    smooth_profile(r->layer_psi6, psi6_s, N_BINS, SMOOTH_HALF);
    smooth_profile(r->layer_rho,  rho_s,  N_BINS, SMOOTH_HALF);

    int gb = -1;
    double p6min = 1e30;
    for(int b = 2; b < N_BINS - 2; b++){
        if(rho_s[b] < RHO_MIN_SMOOTH) continue;
        if(psi6_s[b] < p6min){
            p6min = psi6_s[b];
            gb = b;
        }
    }
    if(gb < 0) gb = N_BINS / 2;

    int bulk = -1;
    double p6max = -1.0;
    for(int b = 2; b < N_BINS - 2; b++){
        if(rho_s[b] < RHO_MIN_SMOOTH) continue;
        int d = abs(b - gb);
        int dwrap = (d < N_BINS - d) ? d : N_BINS - d;
        if(dwrap < MIN_BIN_SEP) continue;
        if(psi6_s[b] > p6max){
            p6max = psi6_s[b];
            bulk = b;
        }
    }
    if(bulk < 0) bulk = (gb + N_BINS / 2) % N_BINS;

    *gb_out   = gb;
    *bulk_out = bulk;
}

static void compute_gb_bulk_averages(MDResult *r, int gb, int bulk)
{
    r->gb_bin   = gb;
    r->bulk_bin = bulk;

    r->psi6_gb    = windowed_average(r->layer_psi6,  r->layer_rho, N_BINS, gb,   GB_WIN_HALF);
    r->psi6_bulk  = windowed_average(r->layer_psi6,  r->layer_rho, N_BINS, bulk, BULK_WIN_HALF);
    r->psi4_gb    = windowed_average(r->layer_psi4,  r->layer_rho, N_BINS, gb,   GB_WIN_HALF);
    r->psi4_bulk  = windowed_average(r->layer_psi4,  r->layer_rho, N_BINS, bulk, BULK_WIN_HALF);
    r->coord_gb   = windowed_average(r->layer_coord, r->layer_rho, N_BINS, gb,   GB_WIN_HALF);
    r->coord_bulk = windowed_average(r->layer_coord, r->layer_rho, N_BINS, bulk, BULK_WIN_HALF);
    r->rho_gb     = windowed_average(r->layer_rho,   r->layer_rho, N_BINS, gb,   GB_WIN_HALF);
    r->rho_bulk   = windowed_average(r->layer_rho,   r->layer_rho, N_BINS, bulk, BULK_WIN_HALF);
}

/* ============================================================ */
static MDResult run_one_temperature(const char *infile, const char *base,
                                    const PointArray *initial,
                                    const TriclinicBox *box_in,
                                    double T_target,
                                    int n_steps, double gamma, double dt,
                                    unsigned long long seed,
                                    int nx, int ny,
                                    int write_outputs,
                                    double P_target, int use_barostat)
{
    MDResult res;
    memset(&res, 0, sizeof(res));
    res.gb_bin   = -1;
    res.bulk_bin = -1;

    int N = initial->n;
    rng_seed(seed);

    PointArray *pos = pa_create(N);
    pa_copy(initial, pos);

    TriclinicBox box = *box_in;

    double *fx = calloc(N, sizeof(double));
    double *fy = calloc(N, sizeof(double));
    double *vx = malloc(N * sizeof(double));
    double *vy = malloc(N * sizeof(double));
    int ncell = nx*ny;
    int *head = malloc(ncell * sizeof(int));
    int *lscl = malloc(N * sizeof(int));

    for(int i=0;i<N;i++){
        vx[i] = sqrt(KB * T_target) * rng_normal();
        vy[i] = sqrt(KB * T_target) * rng_normal();
    }
    double vx_cm = 0.0, vy_cm = 0.0;
    for(int i=0;i<N;i++){ vx_cm += vx[i]; vy_cm += vy[i]; }
    vx_cm /= N; vy_cm /= N;
    for(int i=0;i<N;i++){ vx[i] -= vx_cm; vy[i] -= vy_cm; }

    double virial = 0.0;
    double pe = compute_force_energy(pos, &box, CUTOFF, 1.0, 1.0,
                                     fx, fy, head, lscl, nx, ny, &virial);

    char fname_eng[512], fname_rho[512], fname_prof[512];
    snprintf(fname_eng, sizeof(fname_eng),
             "energy_%s_T%.3f_seed%llu.dat", base, T_target,
             (unsigned long long)seed);
    snprintf(fname_rho, sizeof(fname_rho),
             "rho_%s_T%.3f_seed%llu.dat", base, T_target,
             (unsigned long long)seed);
    snprintf(fname_prof, sizeof(fname_prof),
             "profile_%s_T%.3f_seed%llu.dat", base, T_target,
             (unsigned long long)seed);

    FILE *fp_eng = NULL, *fp_rho = NULL;
    if(write_outputs){
        fp_eng = fopen(fname_eng, "w");
        if(fp_eng){
            fprintf(fp_eng, "# MD run at T = %.6f\n", T_target);
            fprintf(fp_eng, "# input   = %s\n", infile);
            fprintf(fp_eng, "# seed    = %llu\n", (unsigned long long)seed);
            fprintf(fp_eng, "# N_atoms = %d\n", N);
            fprintf(fp_eng, "# barostat= %d  P_target = %.4f\n", use_barostat, P_target);
            fprintf(fp_eng, "# columns : step PE KE T_inst coord psi6 psi4 P\n");
        }
        fp_rho = fopen(fname_rho, "w");
        if(fp_rho){
            fprintf(fp_rho, "# Density profile (accumulated) at T = %.6f\n", T_target);
            fprintf(fp_rho, "# columns : x rho\n");
        }
    }

    double *lrho   = calloc(N_BINS, sizeof(double));
    double *lpsi6  = calloc(N_BINS, sizeof(double));
    double *lpsi4  = calloc(N_BINS, sizeof(double));
    double *lcoord = calloc(N_BINS, sizeof(double));
    int n_layer_samples = 0;

    double coord_sum = 0.0, psi6_sum = 0.0, psi4_sum = 0.0;
    int n_coord_samples = 0;
    double pe_sum = 0.0, ke_sum = 0.0, P_sum = 0.0, a_sum = 0.0;
    int n_samples = 0;

    double c1 = exp(-gamma * dt);
    double c2 = sqrt(KB * T_target * (1.0 - c1*c1));

    int eq_steps = n_steps / 5;
    if(eq_steps < 50000) eq_steps = 50000;
    if(eq_steps > n_steps / 2) eq_steps = n_steps / 2;

    double last_coord = 0.0, last_psi6 = 0.0, last_psi4 = 0.0;

    for(int step=0; step<n_steps; step++){
        for(int i=0;i<N;i++){ vx[i] += 0.5*dt*fx[i]; vy[i] += 0.5*dt*fy[i]; }
        for(int i=0;i<N;i++){ pos->x[i] += 0.5*dt*vx[i]; pos->y[i] += 0.5*dt*vy[i]; }
        for(int i=0;i<N;i++){
            vx[i] = c1*vx[i] + c2*rng_normal();
            vy[i] = c1*vy[i] + c2*rng_normal();
        }
        for(int i=0;i<N;i++){ pos->x[i] += 0.5*dt*vx[i]; pos->y[i] += 0.5*dt*vy[i]; }
        for(int i=0;i<N;i++) wrap_pbc(&box, &pos->x[i], &pos->y[i]);

        pe = compute_force_energy(pos, &box, CUTOFF, 1.0, 1.0,
                                  fx, fy, head, lscl, nx, ny, &virial);

        for(int i=0;i<N;i++){ vx[i] += 0.5*dt*fx[i]; vy[i] += 0.5*dt*fy[i]; }

        if(use_barostat && (step > 0) && (step % BARO_INTERVAL == 0)){
            double P = compute_pressure(&box, vx, vy, N, virial);
            double s = 1.0 + BARO_KAPPA_DEFAULT * (P - P_target) * BARO_INTERVAL * dt;
            if(s < 1.0 - BARO_MAX_SCALE) s = 1.0 - BARO_MAX_SCALE;
            if(s > 1.0 + BARO_MAX_SCALE) s = 1.0 + BARO_MAX_SCALE;
            rescale_box(&box, pos, s);
            pe = compute_force_energy(pos, &box, CUTOFF, 1.0, 1.0,
                                      fx, fy, head, lscl, nx, ny, &virial);
        }

        if(step % 5000 == 0){
            double ke_prog = 0.0;
            for(int i=0;i<N;i++) ke_prog += 0.5*(vx[i]*vx[i] + vy[i]*vy[i]);
            double P_now = compute_pressure(&box, vx, vy, N, virial);
            double a_inst = sqrt(2.0 * box.Lx * box.Ly / (sqrt(3.0) * N));
            printf("      [%s T=%.3f] step %7d/%d  PE/N=%+.4f  T=%.4f  P=%+.4f  a=%.5f\n",
                   base, T_target, step, n_steps, pe/N, ke_prog/N, P_now, a_inst);
            fflush(stdout);
        }

        if(step % SAMPLE_EVERY == 0){
            double ke = 0.0;
            for(int i=0;i<N;i++) ke += 0.5*(vx[i]*vx[i] + vy[i]*vy[i]);
            double T_inst = ke / N;

            double a_inst = sqrt(2.0 * box.Lx * box.Ly / (sqrt(3.0) * N));
            double r_nn = NN_CUTOFF_FRAC * a_inst;

            double coord=0.0, psi6=0.0, psi4=0.0;
            if(step % SAMPLE_COORD == 0){
                double *r_tmp   = calloc(N_BINS, sizeof(double));
                double *p6_tmp  = calloc(N_BINS, sizeof(double));
                double *p4_tmp  = calloc(N_BINS, sizeof(double));
                double *c_tmp   = calloc(N_BINS, sizeof(double));
                accumulate_layer_profile(pos, &box, r_nn,
                                         head, lscl, nx, ny,
                                         r_tmp, p6_tmp, p4_tmp, c_tmp);

                double tot_nb = 0.0;
                for(int b=0;b<N_BINS;b++){
                    lrho[b]   += r_tmp[b];
                    lpsi6[b]  += p6_tmp[b];
                    lpsi4[b]  += p4_tmp[b];
                    lcoord[b] += c_tmp[b];
                    tot_nb    += c_tmp[b];
                }
                n_layer_samples++;

                double p6_w = 0.0, p4_w = 0.0, Nw = 0.0;
                for(int b=0;b<N_BINS;b++){
                    p6_w += p6_tmp[b];
                    p4_w += p4_tmp[b];
                    Nw   += r_tmp[b];
                }
                psi6 = (Nw > 0.0) ? p6_w / Nw : 0.0;
                psi4 = (Nw > 0.0) ? p4_w / Nw : 0.0;
                coord = tot_nb / N;

                last_coord = coord; last_psi6 = psi6; last_psi4 = psi4;
                free(r_tmp); free(p6_tmp); free(p4_tmp); free(c_tmp);
            } else {
                coord = last_coord; psi6 = last_psi6; psi4 = last_psi4;
            }

            double P_inst = compute_pressure(&box, vx, vy, N, virial);

            if(fp_eng){
                fprintf(fp_eng, "%d %.8f %.8f %.6f %.4f %.6f %.6f %.6f\n",
                        step, pe, ke, T_inst, coord, psi6, psi4, P_inst);
            }

            if(step >= eq_steps){
                if(step % SAMPLE_COORD == 0){
                    coord_sum += coord; psi6_sum += psi6; psi4_sum += psi4;
                    n_coord_samples++;
                }
                pe_sum += pe; ke_sum += ke; P_sum += P_inst;
                a_sum  += 0.5*(box.Lx + box.Ly);
                n_samples++;
            }
        }
    }

    if(fp_eng) fclose(fp_eng);

    if(n_layer_samples > 0){
        double bin_w = box.Lx / N_BINS;
        double bin_area = bin_w * box.Ly;
        for(int b=0;b<N_BINS;b++){
            res.layer_rho[b]   = lrho[b] / (n_layer_samples * bin_area);
            double cnt = lrho[b];
            if(cnt > 0.0){
                res.layer_psi6[b]  = lpsi6[b]  / cnt;
                res.layer_psi4[b]  = lpsi4[b]  / cnt;
                res.layer_coord[b] = lcoord[b] / cnt;
            } else {
                res.layer_psi6[b]  = 0.0;
                res.layer_psi4[b]  = 0.0;
                res.layer_coord[b] = 0.0;
            }
        }
        res.n_layer_samples = n_layer_samples;
    }

    res.coord_avg = n_coord_samples ? coord_sum / n_coord_samples : 0.0;
    res.psi6_avg  = n_coord_samples ? psi6_sum  / n_coord_samples : 0.0;
    res.psi4_avg  = n_coord_samples ? psi4_sum  / n_coord_samples : 0.0;
    res.pe_per_atom = n_samples ? pe_sum / (n_samples * N) : 0.0;
    res.e_per_atom  = n_samples ? (pe_sum + ke_sum) / (n_samples * N) : 0.0;
    res.P_avg       = n_samples ? P_sum / n_samples : 0.0;
    res.a_avg       = n_samples ? a_sum / n_samples : 0.0;

    if(write_outputs && fp_rho && n_layer_samples > 0){
        double bin_w = box.Lx / N_BINS;
        for(int b=0;b<N_BINS;b++){
            double x = (b + 0.5) * bin_w;
            fprintf(fp_rho, "%.6f %.8f\n", x, res.layer_rho[b]);
        }
        fclose(fp_rho);
    } else if(fp_rho){ fclose(fp_rho); }

    if(write_outputs && n_layer_samples > 0){
        FILE *fp_prof = fopen(fname_prof, "w");
        if(fp_prof){
            fprintf(fp_prof, "# Layer profile at T = %.6f\n", T_target);
            fprintf(fp_prof, "# input   = %s\n", infile);
            fprintf(fp_prof, "# seed    = %llu\n", (unsigned long long)seed);
            fprintf(fp_prof, "# n_samples = %d\n", n_layer_samples);
            fprintf(fp_prof, "# Lx = %.6f  Ly = %.6f\n", box.Lx, box.Ly);
            fprintf(fp_prof, "# columns : x rho coord psi6 psi4\n");
            double bin_w = box.Lx / N_BINS;
            for(int b=0;b<N_BINS;b++){
                double x = (b + 0.5) * bin_w;
                fprintf(fp_prof, "%.6f %.8f %.6f %.6f %.6f\n",
                        x,
                        res.layer_rho[b],
                        res.layer_coord[b],
                        res.layer_psi6[b],
                        res.layer_psi4[b]);
            }
            fclose(fp_prof);
        }
    }

    free(fx); free(fy); free(vx); free(vy);
    free(head); free(lscl);
    free(lrho); free(lpsi6); free(lpsi4); free(lcoord);
    pa_free(pos);
    return res;
}

/* ============================================================ */
int main(int argc, char **argv){
    printf("=== melt_scan v6.5 (smoothed bin identification, no MSD) ===\n");
    fflush(stdout);

    if(argc < 6){
        fprintf(stderr,
            "Usage: %s <perfect.dat> <bicrystal.dat> <T_min> <T_max> <dT>"
            " [n_steps] [gamma] [dt] [seed] [P_target] [use_barostat]\n"
            "  Fast premelting scan (4 cores):\n"
            "    OMP_NUM_THREADS=4 %s perfect.dat bicrystal.dat "
            "0.28 0.44 0.02 500000 1.0 0.005 5 0.0 1\n",
            argv[0], argv[0]);
        return 1;
    }

    const char *perf_file   = argv[1];
    const char *bicry_file  = argv[2];
    double T_min = atof(argv[3]);
    double T_max = atof(argv[4]);
    double dT    = atof(argv[5]);
    int n_steps      = (argc >  6) ? atoi(argv[6]) : 500000;
    double gamma     = (argc >  7) ? atof(argv[7]) : GAMMA_DEFAULT;
    double dt        = (argc >  8) ? atof(argv[8]) : DT_DEFAULT;
    unsigned long long seed = (argc > 9) ?
        strtoull(argv[9], NULL, 10) : 0xC0FFEEULL;
    double P_target  = (argc > 10) ? atof(argv[10]) : 0.0;
    int use_barostat = (argc > 11) ? atoi(argv[11]) : 1;

    if(dT <= 0.0 || T_max < T_min){
        fprintf(stderr,"Bad T range\n"); return 1;
    }

    int n_threads = omp_get_max_threads();

    printf("[main] reading %s ...\n", perf_file); fflush(stdout);
    PointArray *perf_init = pa_create(1024);
    double p_ax, p_ay, p_bx, p_by, p_a;
    int p_N_hdr;
    if(read_dat(perf_file, perf_init, &p_ax, &p_ay, &p_bx, &p_by, &p_a, &p_N_hdr)){
        pa_free(perf_init); return 1;
    }
    printf("[main]   perfect   N=%d  ax=%.4f by=%.4f  a=%.6f\n",
           perf_init->n, p_ax, p_by, p_a); fflush(stdout);

    TriclinicBox perf_box;
    box_init(&perf_box, p_ax, p_ay, p_bx, p_by,
             -0.5*(p_ax + p_bx), -0.5*(p_ay + p_by));

    printf("[main] reading %s ...\n", bicry_file); fflush(stdout);
    PointArray *bicry_init = pa_create(1024);
    double b_ax, b_ay, b_bx, b_by, b_a;
    int b_N_hdr;
    if(read_dat(bicry_file, bicry_init, &b_ax, &b_ay, &b_bx, &b_by, &b_a, &b_N_hdr)){
        pa_free(perf_init); pa_free(bicry_init); return 1;
    }
    printf("[main]   bicrystal N=%d  ax=%.4f by=%.4f  a=%.6f\n",
           bicry_init->n, b_ax, b_by, b_a); fflush(stdout);

    TriclinicBox bicry_box;
    box_init(&bicry_box, b_ax, b_ay, b_bx, b_by,
             -0.5*(b_ax + b_bx), -0.5*(b_ay + b_by));

    double phA = fabs(perf_box.det) / perf_box.Ly;
    double phB = fabs(perf_box.det) / perf_box.Lx;
    int p_nx = (int)(phA / CUTOFF); if(p_nx < 1) p_nx = 1;
    int p_ny = (int)(phB / CUTOFF); if(p_ny < 1) p_ny = 1;

    double bhA = fabs(bicry_box.det) / bicry_box.Ly;
    double bhB = fabs(bicry_box.det) / bicry_box.Lx;
    int b_nx = (int)(bhA / CUTOFF); if(b_nx < 1) b_nx = 1;
    int b_ny = (int)(bhB / CUTOFF); if(b_ny < 1) b_ny = 1;

    char p_base[128], b_base[128];
    basename_no_ext(perf_file,  p_base, sizeof(p_base));
    basename_no_ext(bicry_file, b_base, sizeof(b_base));

    int n_temps = (int)floor((T_max - T_min) / dT + 1e-9) + 1;

    printf("\n==========================================================\n");
    printf("  T range    : %.4f .. %.4f  (dT = %.4f, %d points)\n",
           T_min, T_max, dT, n_temps);
    printf("  n_steps    : %d\n", n_steps);
    printf("  gamma / dt : %.3f / %.4f\n", gamma, dt);
    printf("  seed       : %llu\n", (unsigned long long)seed);
    printf("  P_target   : %.4f  (barostat %s)\n",
           P_target, use_barostat ? "ON" : "OFF");
    printf("  OMP threads: %d\n", n_threads);
    printf("==========================================================\n\n");
    fflush(stdout);

    {
        int gb0 = -1, bulk0 = -1;
        identify_gb_and_bulk_bins(bicry_init, &bicry_box,
                                  b_nx, b_ny,
                                  &gb0, &bulk0);
        printf("\n");
    }

    MDResult *results_p = calloc(n_temps, sizeof(MDResult));
    MDResult *results_b = calloc(n_temps, sizeof(MDResult));

    #pragma omp parallel for schedule(dynamic)
    for(int k=0; k<n_temps; k++){
        double T = T_min + k * dT;
        unsigned long long seed_T = seed + (unsigned long long)(k * 1000003ULL);

        results_p[k] = run_one_temperature(perf_file, p_base,
                                           perf_init, &perf_box,
                                           T, n_steps, gamma, dt,
                                           seed_T + 0ULL,
                                           p_nx, p_ny, 1,
                                           P_target, use_barostat);
        results_b[k] = run_one_temperature(bicry_file, b_base,
                                           bicry_init, &bicry_box,
                                           T, n_steps, gamma, dt,
                                           seed_T + 777ULL,
                                           b_nx, b_ny, 1,
                                           P_target, use_barostat);
    }

    printf("\n==========================================================\n");
    printf("  Post-processing: identifying GB/bulk bins from smoothed profiles\n");
    printf("==========================================================\n");
    for(int k=0; k<n_temps; k++){
        int gb, bulk;
        identify_bins_from_profile(&results_b[k], &gb, &bulk);
        compute_gb_bulk_averages(&results_b[k], gb, bulk);
        compute_gb_bulk_averages(&results_p[k], gb, bulk);

        printf("  T=%.3f  gb_bin=%2d  bulk_bin=%2d  "
               "(GB psi6_B=%.4f, bulk psi6_B=%.4f, bulk psi6_P=%.4f)\n",
               T_min + k * dT, gb, bulk,
               results_b[k].psi6_gb, results_b[k].psi6_bulk,
               results_p[k].psi6_bulk);
    }
    fflush(stdout);
    printf("\n");

    char fname_sum[512];
    snprintf(fname_sum, sizeof(fname_sum),
             "sweep_%s_vs_%s.dat", p_base, b_base);
    FILE *fp_sum = fopen(fname_sum, "w");
    if(fp_sum){
        fprintf(fp_sum, "# T-scan comparison (v6.5, smoothed bin id)\n");
        fprintf(fp_sum, "# perfect   = %s\n", perf_file);
        fprintf(fp_sum, "# bicrystal = %s\n", bicry_file);
        fprintf(fp_sum, "# n_steps   = %d\n", n_steps);
        fprintf(fp_sum, "# gamma     = %.6f\n", gamma);
        fprintf(fp_sum, "# dt        = %.6f\n", dt);
        fprintf(fp_sum, "# seed      = %llu\n", (unsigned long long)seed);
        fprintf(fp_sum, "# P_target  = %.6f\n", P_target);
        fprintf(fp_sum, "# barostat  = %d\n", use_barostat);
        fprintf(fp_sum, "# omp_threads = %d\n", n_threads);
        fprintf(fp_sum, "# gb_win_half    = %d\n", GB_WIN_HALF);
        fprintf(fp_sum, "# bulk_win_half  = %d\n", BULK_WIN_HALF);
        fprintf(fp_sum, "# min_bin_sep    = %d\n", MIN_BIN_SEP);
        fprintf(fp_sum, "# smooth_half    = %d\n", SMOOTH_HALF);
        fprintf(fp_sum, "# NOTE: gb_bin and bulk_bin are identified per-T from\n");
        fprintf(fp_sum, "#       the smoothed bicrystal production-averaged profile.\n");
        fprintf(fp_sum,
            "# columns : T "
            "psi6_avg_P psi6_gb_P psi6_bulk_P PE_P P_P "
            "psi6_avg_B psi6_gb_B psi6_bulk_B PE_B P_B "
            "gb_bin_B bulk_bin_B "
            "d_psi6_avg d_psi6_gb d_psi6_bulk d_PE\n");
    }

    for(int k=0; k<n_temps; k++){
        double T = T_min + k * dT;
        MDResult rp = results_p[k];
        MDResult rb = results_b[k];

        printf("----------------------------------------------------------\n");
        printf("  T = %.4f  (gb_bin=%d, bulk_bin=%d)\n",
               T, rb.gb_bin, rb.bulk_bin);
        printf("  perfect   : psi6 avg=%.4f gb=%.4f bulk=%.4f  PE/N=%.6f  P=%.4f\n",
               rp.psi6_avg, rp.psi6_gb, rp.psi6_bulk,
               rp.pe_per_atom, rp.P_avg);
        printf("  bicrystal : psi6 avg=%.4f gb=%.4f bulk=%.4f  PE/N=%.6f  P=%.4f\n",
               rb.psi6_avg, rb.psi6_gb, rb.psi6_bulk,
               rb.pe_per_atom, rb.P_avg);
        fflush(stdout);

        double d_psi6_avg  = rb.psi6_avg  - rp.psi6_avg;
        double d_psi6_gb   = rb.psi6_gb   - rp.psi6_gb;
        double d_psi6_bulk = rb.psi6_bulk - rp.psi6_bulk;
        double d_PE        = rb.pe_per_atom - rp.pe_per_atom;

        if(fp_sum){
            fprintf(fp_sum,
                "%.4f "
                "%.6f %.6f %.6f %.6f %.6f "
                "%.6f %.6f %.6f %.6f %.6f "
                "%d %d "
                "%+.6f %+.6f %+.6f %+.6f\n",
                T,
                rp.psi6_avg, rp.psi6_gb, rp.psi6_bulk, rp.pe_per_atom, rp.P_avg,
                rb.psi6_avg, rb.psi6_gb, rb.psi6_bulk, rb.pe_per_atom, rb.P_avg,
                rb.gb_bin, rb.bulk_bin,
                d_psi6_avg, d_psi6_gb, d_psi6_bulk, d_PE);
        }
    }

    if(fp_sum) fclose(fp_sum);

    printf("\n==========================================================\n");
    printf("Sweep summary -> %s\n", fname_sum);
    printf("==========================================================\n");

    free(results_p);
    free(results_b);
    pa_free(perf_init);
    pa_free(bicry_init);
    return 0;
}