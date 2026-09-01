#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

/* ============================================================
   PARAMETERS
   ============================================================ */
#define DT          0.005
#define CUTOFF      3.2
#define CUTOFF2     (CUTOFF * CUTOFF)
#define KB          1.0

/* Box for 100a x 30a with a = 1.2 */
#define BOX_XMIN    -120.0
#define BOX_XMAX     120.0
#define BOX_YMIN     -36.373
#define BOX_YMAX      36.373
#define LX          (BOX_XMAX - BOX_XMIN)
#define LY          (BOX_YMAX - BOX_YMIN)

#define T_START     0.05
#define T_END       1.0
#define DTEMP       0.01
#define N_EQ        1000
#define N_PROD      2000

#define TAU         0.1
#define PI          3.141592653589793
#define NEIGHBOR_CUT 1.8
#define NEIGHBOR_CUT2 (NEIGHBOR_CUT * NEIGHBOR_CUT)

#define N_SAMPLES   10

/* GB probe region: |x| < GB_LIMIT */
#define GB_LIMIT    3.0

/* Fixed random seed used in the code */
#define SEED        6

/* ============================================================
   Precomputed shifted‑LJ cutoff constants
   ============================================================ */
static const double RC           = CUTOFF;
static const double RC_INV       = 1.0 / CUTOFF;
static const double RC2_INV      = RC_INV * RC_INV;
static const double RC6_INV      = RC2_INV * RC2_INV * RC2_INV;
static const double RC12_INV     = RC6_INV * RC6_INV;
static const double U_RC         = 4.0 * (RC12_INV - RC6_INV);
static const double DUDR_RC      = -24.0 * (2.0 * RC12_INV - RC6_INV) * RC_INV;
static const double F_RC         = -DUDR_RC;
static const double SHIFT        = U_RC;

/* ============================================================
   Structure: Atoms
   ============================================================ */
typedef struct {
    int N;
    double *x, *y;
    double *vx, *vy;
    double *fx, *fy;
} Atoms;

/* ============================================================
   Shifted‑force Lennard‑Jones potential and force
   ============================================================ */
static inline void lj_shifted_force_pe(double dx, double dy,
                                       double *fx, double *fy, double *pe) {
    double r2 = dx*dx + dy*dy;
    if (r2 >= CUTOFF2 || r2 < 1e-12) {
        *fx = 0.0; *fy = 0.0; *pe = 0.0;
        return;
    }

    double r = sqrt(r2);
    double r2inv = 1.0 / r2;
    double r6inv = r2inv * r2inv * r2inv;
    double r12inv = r6inv * r6inv;

    double f_raw = 24.0 * (2.0 * r12inv - r6inv) / r;
    double f_eff = f_raw - F_RC;
    double f = f_eff / r;

    double U = 4.0 * (r12inv - r6inv);
    *fx = f * dx;
    *fy = f * dy;
    *pe = U - SHIFT - F_RC * (RC - r);
}

/* ============================================================
   Compute forces and potential energy — CELL LIST VERSION
   ============================================================ */
void compute_force_and_pe(Atoms *atoms, double *pe_out) {
    int N = atoms->N;
    double *x = atoms->x, *y = atoms->y;
    double *fx = atoms->fx, *fy = atoms->fy;

    for (int i = 0; i < N; i++) {
        fx[i] = 0.0;
        fy[i] = 0.0;
    }

    double pe = 0.0;

    double cell_size = CUTOFF;
    int nx = (int)ceil(LX / cell_size);
    int ny = (int)ceil(LY / cell_size);
    if (nx < 1) nx = 1;
    if (ny < 1) ny = 1;
    int ncell = nx * ny;

    int *head = malloc(ncell * sizeof(int));
    int *lscl = malloc(N * sizeof(int));
    if (!head || !lscl) {
        fprintf(stderr, "Allocation error in cell list\n");
        exit(1);
    }
    for (int c = 0; c < ncell; c++) head[c] = -1;

    for (int i = 0; i < N; i++) {
        int cx = (int)floor((x[i] - BOX_XMIN) / cell_size);
        if (cx < 0) cx = 0;
        if (cx >= nx) cx = nx - 1;

        int cy = (int)floor((y[i] - BOX_YMIN) / cell_size);
        cy = cy % ny;
        if (cy < 0) cy += ny;

        int cell = cy * nx + cx;
        lscl[i] = head[cell];
        head[cell] = i;
    }

    double dx, dy, fx_pair, fy_pair, pe_pair;
    for (int i = 0; i < N; i++) {
        int cx = (int)floor((x[i] - BOX_XMIN) / cell_size);
        if (cx < 0) cx = 0;
        if (cx >= nx) cx = nx - 1;
        int cy = (int)floor((y[i] - BOX_YMIN) / cell_size);
        cy = cy % ny;
        if (cy < 0) cy += ny;

        for (int ii = cx - 1; ii <= cx + 1; ii++) {
            if (ii < 0 || ii >= nx) continue;
            for (int jj = cy - 1; jj <= cy + 1; jj++) {
                int ccy = jj % ny;
                if (ccy < 0) ccy += ny;
                int cell = ccy * nx + ii;

                for (int j = head[cell]; j != -1; j = lscl[j]) {
                    if (j <= i) continue;

                    dx = x[i] - x[j];
                    dy = y[i] - y[j];
                    dy -= LY * round(dy / LY);

                    lj_shifted_force_pe(dx, dy, &fx_pair, &fy_pair, &pe_pair);

                    fx[i] += fx_pair;
                    fy[i] += fy_pair;
                    fx[j] -= fx_pair;
                    fy[j] -= fy_pair;
                    pe += pe_pair;
                }
            }
        }
    }

    free(head);
    free(lscl);
    *pe_out = pe;
}

/* ============================================================
   Velocity Verlet + PBC in y + REFLECTING WALLS IN X
   ============================================================ */
void velocity_verlet(Atoms *atoms) {
    int N = atoms->N;
    double *x = atoms->x, *y = atoms->y;
    double *vx = atoms->vx, *vy = atoms->vy;
    double *fx = atoms->fx, *fy = atoms->fy;

    for (int i = 0; i < N; i++) {
        vx[i] += 0.5 * DT * fx[i];
        vy[i] += 0.5 * DT * fy[i];
    }

    for (int i = 0; i < N; i++) {
        x[i] += vx[i] * DT;
        y[i] += vy[i] * DT;
    }

    // PBC in y
    for (int i = 0; i < N; i++) {
        y[i] = BOX_YMIN + fmod(y[i] - BOX_YMIN, LY);
        if (y[i] < BOX_YMIN) y[i] += LY;
    }

    // Reflecting walls in x
    for (int i = 0; i < N; i++) {
        if (x[i] < BOX_XMIN) {
            x[i] = 2.0 * BOX_XMIN - x[i];
            vx[i] = -vx[i];
        } else if (x[i] > BOX_XMAX) {
            x[i] = 2.0 * BOX_XMAX - x[i];
            vx[i] = -vx[i];
        }
    }

    double pe;
    compute_force_and_pe(atoms, &pe);

    for (int i = 0; i < N; i++) {
        vx[i] += 0.5 * DT * atoms->fx[i];
        vy[i] += 0.5 * DT * atoms->fy[i];
    }
}

/* ============================================================
   Berendsen Thermostat (correct 2D temperature)
   ============================================================ */
void berendsen_thermostat(Atoms *atoms, double target_T) {
    int N = atoms->N;
    double *vx = atoms->vx, *vy = atoms->vy;
    double KE = 0.0;
    for (int i = 0; i < N; i++) {
        KE += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i]);
    }
    double current_T = KE / (N * KB);
    if (current_T > 0.0) {
        double factor = 1.0 + (DT / TAU) * (target_T / current_T - 1.0);
        if (factor < 0.0) factor = 0.0;
        double lambda = sqrt(factor);
        for (int i = 0; i < N; i++) {
            vx[i] *= lambda;
            vy[i] *= lambda;
        }
    }
}

/* ============================================================
   Compute local |ψ6| for all atoms — CELL LIST + OpenMP
   ============================================================ */
double *compute_psi6_all(Atoms *atoms) {
    int N = atoms->N;
    double *x = atoms->x, *y = atoms->y;
    double *psi6 = malloc(N * sizeof(double));
    if (!psi6) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    double cell_size = NEIGHBOR_CUT;
    int nx = (int)ceil(LX / cell_size);
    int ny = (int)ceil(LY / cell_size);
    if (nx < 1) nx = 1;
    if (ny < 1) ny = 1;
    int ncell = nx * ny;

    int *head = malloc(ncell * sizeof(int));
    int *lscl = malloc(N * sizeof(int));
    if (!head || !lscl) {
        fprintf(stderr, "Allocation error in psi6 cell list\n");
        exit(1);
    }
    for (int c = 0; c < ncell; c++) head[c] = -1;

    for (int i = 0; i < N; i++) {
        int cx = (int)floor((x[i] - BOX_XMIN) / cell_size);
        if (cx < 0) cx = 0;
        if (cx >= nx) cx = nx - 1;

        int cy = (int)floor((y[i] - BOX_YMIN) / cell_size);
        cy = cy % ny;
        if (cy < 0) cy += ny;

        int cell = cy * nx + cx;
        lscl[i] = head[cell];
        head[cell] = i;
    }

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int i = 0; i < N; i++) {
        double real = 0.0, imag = 0.0;
        int neigh = 0;

        int cx = (int)floor((x[i] - BOX_XMIN) / cell_size);
        if (cx < 0) cx = 0;
        if (cx >= nx) cx = nx - 1;
        int cy = (int)floor((y[i] - BOX_YMIN) / cell_size);
        cy = cy % ny;
        if (cy < 0) cy += ny;

        for (int ii = cx - 1; ii <= cx + 1; ii++) {
            if (ii < 0 || ii >= nx) continue;
            for (int jj = cy - 1; jj <= cy + 1; jj++) {
                int ccy = jj % ny;
                if (ccy < 0) ccy += ny;
                int cell = ccy * nx + ii;

                for (int j = head[cell]; j != -1; j = lscl[j]) {
                    if (j == i) continue;

                    double dx = x[i] - x[j];
                    double dy = y[i] - y[j];
                    dy -= LY * round(dy / LY);
                    double r2 = dx*dx + dy*dy;

                    if (r2 < NEIGHBOR_CUT2 && r2 > 1e-12) {
                        double theta = atan2(dy, dx);
                        real += cos(6.0 * theta);
                        imag += sin(6.0 * theta);
                        neigh++;
                    }
                }
            }
        }

        if (neigh > 0) {
            real /= neigh;
            imag /= neigh;
            psi6[i] = sqrt(real*real + imag*imag);
        } else {
            psi6[i] = 0.0;
        }
    }

    free(head);
    free(lscl);
    return psi6;
}

/* ============================================================
   Initialize velocities from Maxwell–Boltzmann at temperature T
   ============================================================ */
void initialize_velocities(Atoms *atoms, double T) {
    int N = atoms->N;
    for (int i = 0; i < N; i++) {
        double u1 = (rand() + 1.0) / (RAND_MAX + 2.0);
        double u2 = (rand() + 1.0) / (RAND_MAX + 2.0);
        double mag = sqrt(-2.0 * T * log(u1));
        double angle = 2.0 * PI * u2;
        atoms->vx[i] = mag * cos(angle);
        atoms->vy[i] = mag * sin(angle);
    }

    // Remove net momentum
    double sum_vx = 0.0, sum_vy = 0.0;
    for (int i = 0; i < N; i++) {
        sum_vx += atoms->vx[i];
        sum_vy += atoms->vy[i];
    }
    sum_vx /= N;
    sum_vy /= N;
    for (int i = 0; i < N; i++) {
        atoms->vx[i] -= sum_vx;
        atoms->vy[i] -= sum_vy;
    }

    // Rescale to exactly match target T after momentum removal
    double KE = 0.0;
    for (int i = 0; i < N; i++) {
        KE += 0.5 * (atoms->vx[i]*atoms->vx[i] + atoms->vy[i]*atoms->vy[i]);
    }
    double current_T = KE / (N * KB);
    if (current_T > 0.0) {
        double scale = sqrt(T / current_T);
        for (int i = 0; i < N; i++) {
            atoms->vx[i] *= scale;
            atoms->vy[i] *= scale;
        }
    }
}

/* ============================================================
   Read atoms from file — ignores extra columns, with safety checks
   ============================================================ */
Atoms read_atoms(const char *filename) {
    Atoms atoms;
    atoms.N = 0;
    atoms.x = NULL;
    atoms.y = NULL;
    atoms.vx = NULL;
    atoms.vy = NULL;
    atoms.fx = NULL;
    atoms.fy = NULL;

    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("Error opening file");
        exit(1);
    }

    double *raw_x = NULL, *raw_y = NULL;
    int N_total = 0;
    char line[512];
    double tx, ty;

    while (fgets(line, sizeof(line), fp)) {
        if (sscanf(line, "%lf %lf", &tx, &ty) == 2) {
            double *new_x = realloc(raw_x, (N_total + 1) * sizeof(double));
            if (!new_x) {
                free(raw_x);
                free(raw_y);
                fclose(fp);
                fprintf(stderr, "realloc failed while reading %s\n", filename);
                exit(1);
            }
            raw_x = new_x;

            double *new_y = realloc(raw_y, (N_total + 1) * sizeof(double));
            if (!new_y) {
                free(raw_x);
                free(raw_y);
                fclose(fp);
                fprintf(stderr, "realloc failed while reading %s\n", filename);
                exit(1);
            }
            raw_y = new_y;

            raw_x[N_total] = tx;
            raw_y[N_total] = ty;
            N_total++;
        }
    }
    fclose(fp);

    if (N_total == 0) {
        fprintf(stderr, "No atoms found in %s\n", filename);
        exit(1);
    }

    atoms.N = N_total;
    atoms.x = malloc(N_total * sizeof(double));
    atoms.y = malloc(N_total * sizeof(double));
    atoms.vx = calloc(N_total, sizeof(double));
    atoms.vy = calloc(N_total, sizeof(double));
    atoms.fx = malloc(N_total * sizeof(double));
    atoms.fy = malloc(N_total * sizeof(double));

    if (!atoms.x || !atoms.y || !atoms.vx || !atoms.vy || !atoms.fx || !atoms.fy) {
        fprintf(stderr, "Memory allocation failed for atoms structure\n");
        exit(1);
    }

    for (int i = 0; i < N_total; i++) {
        atoms.x[i] = raw_x[i];
        atoms.y[i] = raw_y[i];
        atoms.y[i] = BOX_YMIN + fmod(atoms.y[i] - BOX_YMIN, LY);
        if (atoms.y[i] < BOX_YMIN) atoms.y[i] += LY;
    }

    free(raw_x);
    free(raw_y);
    return atoms;
}

/* ============================================================
   Free atoms structure
   ============================================================ */
void free_atoms(Atoms *atoms) {
    free(atoms->x);
    free(atoms->y);
    free(atoms->vx);
    free(atoms->vy);
    free(atoms->fx);
    free(atoms->fy);
}

/* ============================================================
   Run simulation, return average GB ψ6 for each temperature
   ============================================================ */
double *run_simulation(const char *filename, double **T_values, int *nT) {
    Atoms atoms = read_atoms(filename);
    initialize_velocities(&atoms, T_START);

    double pe;
    compute_force_and_pe(&atoms, &pe);

    int count_T = 0;
    for (double T = T_START; T <= T_END + 1e-6; T += DTEMP) count_T++;

    double *gb_avg = malloc(count_T * sizeof(double));
    double *temp_arr = malloc(count_T * sizeof(double));
    if (!gb_avg || !temp_arr) {
        fprintf(stderr, "Memory allocation failed in run_simulation\n");
        free_atoms(&atoms);
        exit(1);
    }

    int T_idx = 0;
    double T = T_START;

    while (T <= T_END + 1e-6) {
        printf("  T = %.2f: Equilibrating...\n", T);
        for (int step = 0; step < N_EQ; step++) {
            velocity_verlet(&atoms);
            berendsen_thermostat(&atoms, T);
        }

        printf("  T = %.2f: Production...\n", T);

        int sample_interval = N_PROD / N_SAMPLES;
        int sample_idx = 0;
        int step_counter = 0;
        double psi_sum = 0.0;
        int psi_count = 0;

        for (int step = 0; step < N_PROD; step++) {
            velocity_verlet(&atoms);
            berendsen_thermostat(&atoms, T);

            step_counter++;
            if (step_counter >= sample_interval && sample_idx < N_SAMPLES) {
                double *psi6 = compute_psi6_all(&atoms);

                double gb_sum = 0.0;
                int gb_count = 0;
                for (int i = 0; i < atoms.N; i++) {
                    if (fabs(atoms.x[i]) < GB_LIMIT) {
                        gb_sum += psi6[i];
                        gb_count++;
                    }
                }
                if (gb_count > 0) {
                    psi_sum += gb_sum / gb_count;
                    psi_count++;
                } else {
                    printf("    Warning: no atoms in GB probe at sample %d\n", sample_idx);
                }

                free(psi6);
                sample_idx++;
                step_counter = 0;
            }
        }

        if (psi_count > 0) {
            gb_avg[T_idx] = psi_sum / psi_count;
        } else {
            if (T_idx > 0) {
                gb_avg[T_idx] = gb_avg[T_idx - 1] - 0.001;
                printf("  Warning: no valid GB sample at T=%.2f, using previous value.\n", T);
            } else {
                gb_avg[T_idx] = 0.0;
            }
        }

        printf("  T = %.2f: GB |ψ6| = %.6f\n", T, gb_avg[T_idx]);

        temp_arr[T_idx] = T;
        T_idx++;
        T += DTEMP;
    }

    free_atoms(&atoms);

    *T_values = temp_arr;
    *nT = count_T;
    return gb_avg;
}

/* ============================================================
   Helper: make output filename from input file basename
   ============================================================ */
void make_output_filename(const char *input_path, char *out, size_t out_size) {
    const char *base = input_path;

    // Strip directory (both Unix and Windows separators)
    const char *slash = strrchr(input_path, '/');
    const char *backslash = strrchr(input_path, '\\');
    if (slash && slash > base) base = slash + 1;
    if (backslash && backslash > base) base = backslash + 1;

    // Copy stem
    char stem[256];
    strncpy(stem, base, sizeof(stem) - 1);
    stem[sizeof(stem) - 1] = '\0';

    // Remove extension
    char *dot = strrchr(stem, '.');
    if (dot) *dot = '\0';

    snprintf(out, out_size, "%s_seed%d_psi_GB_vs_T.csv", stem, SEED);
}

/* ============================================================
   Main — multiple input files version
   ============================================================ */
int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <input_file1> [input_file2] ...\n", argv[0]);
        return 1;
    }

    // Fixed seed for reproducibility
    srand(SEED);

    clock_t global_start = clock();

    for (int f = 1; f < argc; f++) {
        printf("\n============================================================\n");
        printf("Running simulation for: %s\n", argv[f]);

        double *T_arr = NULL;
        int nT = 0;
        double *gb_psi = run_simulation(argv[f], &T_arr, &nT);

        char out_filename[512];
        make_output_filename(argv[f], out_filename, sizeof(out_filename));

        FILE *out = fopen(out_filename, "w");
        if (!out) {
            perror("Error opening output file");
            free(gb_psi);
            free(T_arr);
            continue;
        }
        fprintf(out, "T,psi_GB\n");

        printf("\n=== FINAL GB |ψ6| TABLE for %s ===\n", argv[f]);
        printf("T        psi_GB\n");
        for (int t = 0; t < nT; t++) {
            printf("%.2f    %.6f\n", T_arr[t], gb_psi[t]);
            fprintf(out, "%.6f,%.6f\n", T_arr[t], gb_psi[t]);
        }

        fclose(out);
        free(gb_psi);
        free(T_arr);

        printf("GB |ψ6| CSV written to: %s\n", out_filename);
    }

    clock_t global_end = clock();
    double elapsed = (double)(global_end - global_start) / CLOCKS_PER_SEC;
    printf("\nAll simulations finished in %.2f seconds (%.2f minutes).\n", elapsed, elapsed/60.0);

    return 0;
}