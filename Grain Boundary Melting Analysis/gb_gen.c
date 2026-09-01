#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <omp.h>

#define PI 3.14159265358979323846
#define EPS 1e-12

/* ---------- dynamic array of 2D points ---------- */
typedef struct {
    double *x, *y;
    int n, cap;
} PointArray;

PointArray* pa_create(int cap) {
    PointArray *a = malloc(sizeof(PointArray));
    if (!a) return NULL;
    a->x = malloc(cap * sizeof(double));
    a->y = malloc(cap * sizeof(double));
    a->n = 0;
    a->cap = cap;
    return a;
}

void pa_free(PointArray *a) {
    if (a) { free(a->x); free(a->y); free(a); }
}

void pa_append(PointArray *a, double x, double y) {
    if (a->n >= a->cap) {
        a->cap = a->cap ? a->cap * 2 : 8;
        a->x = realloc(a->x, a->cap * sizeof(double));
        a->y = realloc(a->y, a->cap * sizeof(double));
    }
    a->x[a->n] = x;
    a->y[a->n] = y;
    a->n++;
}

void pa_copy(PointArray *src, PointArray *dst) {
    dst->n = 0;
    for (int i = 0; i < src->n; i++)
        pa_append(dst, src->x[i], src->y[i]);
}

void pa_clear(PointArray *a) { a->n = 0; }

/* ---------- exact periodic y wrap ---------- */
static inline double wrap_y(double y, double ymin, double Ly) {
    double t = fmod(y - ymin, Ly);
    if (t < 0.0)
        t += Ly;
    return ymin + t;
}

/* ---------- lattice generation ---------- */
PointArray* generate_triangular_lattice(double a,
                                        double xmin, double xmax,
                                        double ymin, double ymax,
                                        double margin) {
    PointArray *pts = pa_create(1000);
    double e1x = a, e1y = 0.0;
    double e2x = a * 0.5, e2y = a * sqrt(3.0) / 2.0;
    int N = (int)(fmax(xmax - xmin + 2*margin, ymax - ymin + 2*margin) / a) + 5;
    for (int i = -N; i <= N; i++) {
        for (int j = -N; j <= N; j++) {
            double x = i * e1x + j * e2x;
            double y = i * e1y + j * e2y;
            if (x >= xmin - margin && x <= xmax + margin &&
                y >= ymin - margin && y <= ymax + margin)
                pa_append(pts, x, y);
        }
    }
    return pts;
}

void filter_box(PointArray *src, double xmin, double xmax,
                double ymin, double ymax, PointArray *dst) {
    pa_clear(dst);
    for (int i = 0; i < src->n; i++) {
        double x = src->x[i], y = src->y[i];
        if (x >= xmin && x <= xmax && y >= ymin && y <= ymax)
            pa_append(dst, x, y);
    }
}

void rotate_points(PointArray *src, double theta, PointArray *dst) {
    double c = cos(theta), s = sin(theta);
    pa_clear(dst);
    for (int i = 0; i < src->n; i++) {
        double x = src->x[i], y = src->y[i];
        pa_append(dst, c*x - s*y, s*x + c*y);
    }
}

/* ============================================================
   CELL‑LIST BASED LJ FORCE/ENERGY (O(N) instead of O(N²))
   ============================================================ */
double lj_force_energy_cell(PointArray *pos,
                            double box_xmin, double box_xmax,
                            double box_ymin, double box_ymax,
                            double rc, double sigma, double epsilon,
                            double *fx, double *fy) {
    int N = pos->n;
    double rc2 = rc * rc;
    double Lx = box_xmax - box_xmin;
    double Ly = box_ymax - box_ymin;
    double energy = 0.0;

    double rc6 = pow(sigma / rc, 6.0);
    double rc12 = rc6 * rc6;
    double shift = 4.0 * epsilon * (rc12 - rc6);
    double f_rc = 4.0 * epsilon * (12.0 * rc12 - 6.0 * rc6) / rc;

    memset(fx, 0, N * sizeof(double));
    memset(fy, 0, N * sizeof(double));

    // Cell list setup
    double cell_size = rc;               // must be >= rc
    int nx = (int) ceil(Lx / cell_size);
    int ny = (int) ceil(Ly / cell_size);
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

    // Build linked cell list
    for (int i = 0; i < N; i++) {
        double xi = pos->x[i];
        double yi = pos->y[i];

        int cx = (int) floor((xi - box_xmin) / cell_size);
        if (cx < 0) cx = 0;
        if (cx >= nx) cx = nx - 1;

        int cy = (int) floor((yi - box_ymin) / cell_size);
        cy = cy % ny;
        if (cy < 0) cy += ny;

        int cell = cy * nx + cx;
        lscl[i] = head[cell];
        head[cell] = i;
    }

    // Loop over atoms and neighbouring cells
    for (int i = 0; i < N; i++) {
        double xi = pos->x[i];
        double yi = pos->y[i];

        int cx = (int) floor((xi - box_xmin) / cell_size);
        if (cx < 0) cx = 0;
        if (cx >= nx) cx = nx - 1;
        int cy = (int) floor((yi - box_ymin) / cell_size);
        cy = cy % ny;
        if (cy < 0) cy += ny;

        for (int ii = cx - 1; ii <= cx + 1; ii++) {
            int ccx = ii;
            if (ccx < 0 || ccx >= nx) continue;   // no PBC in x
            for (int jj = cy - 1; jj <= cy + 1; jj++) {
                int ccy = jj % ny;
                if (ccy < 0) ccy += ny;
                int cell = ccy * nx + ccx;

                for (int j = head[cell]; j != -1; j = lscl[j]) {
                    if (j <= i) continue;   // avoid double counting & self

                    double dx = xi - pos->x[j];
                    double dy = yi - pos->y[j];
                    dy -= Ly * round(dy / Ly);   // PBC in y

                    double r2 = dx*dx + dy*dy;
                    if (r2 < rc2 && r2 > 1e-12) {
                        double r = sqrt(r2);
                        double r6 = pow(sigma / r, 6.0);
                        double r12 = r6 * r6;

                        double f_raw = 4.0 * epsilon * (12.0 * r12 - 6.0 * r6) / r;
                        double f_phys_eff = f_raw - f_rc;
                        double f_mag = f_phys_eff / r;

                        double fxi = f_mag * dx;
                        double fyi = f_mag * dy;
                        fx[i] += fxi;  fy[i] += fyi;
                        fx[j] -= fxi;  fy[j] -= fyi;

                        energy += 4.0 * epsilon * (r12 - r6)
                                  - shift - f_rc * (rc - r);
                    }
                }
            }
        }
    }

    free(head);
    free(lscl);
    return energy;
}

/* ---------- gradient descent relaxation with displacement cap ---------- */
void gradient_descent_relax(PointArray *pos,
                            double box_xmin, double box_xmax,
                            double box_ymin, double box_ymax,
                            double rc, double sigma, double epsilon,
                            double dt, double damping, int max_steps, double ftol,
                            bool verbose,
                            bool *fixed_mask, PointArray *fixed_positions,
                            double max_disp) {
    int N = pos->n;
    double *fx = malloc(N * sizeof(double));
    double *fy = malloc(N * sizeof(double));
    double Ly = box_ymax - box_ymin;

    PointArray *fixed_copy = pa_create(N);
    if (fixed_positions) pa_copy(fixed_positions, fixed_copy);

    for (int step = 0; step < max_steps; step++) {
        double energy = lj_force_energy_cell(pos, box_xmin, box_xmax,
                                             box_ymin, box_ymax,
                                             rc, sigma, epsilon, fx, fy);
        double maxf = 0.0;
        for (int i = 0; i < N; i++) {
            double f = sqrt(fx[i]*fx[i] + fy[i]*fy[i]);
            if (f > maxf) maxf = f;
        }

        if (verbose && step % 100 == 0)
            printf("  Step %5d, max force = %.2e, energy = %.6f\n",
                   step, maxf, energy);

        if (maxf < ftol) {
            if (verbose) printf("  Converged after %d steps.\n", step);
            break;
        }

        // Gradient descent step
        for (int i = 0; i < N; i++) {
            if (fixed_mask && fixed_mask[i]) {
                pos->x[i] = fixed_copy->x[i];
                pos->y[i] = fixed_copy->y[i];
                continue;
            }

            double dx = dt * fx[i];
            double dy = dt * fy[i];

            double d = sqrt(dx*dx + dy*dy);
            if (d > max_disp) {
                double scale = max_disp / d;
                dx *= scale;
                dy *= scale;
            }

            pos->x[i] += dx;
            pos->y[i] += dy;
        }

        // Periodic wrap in y
        for (int i = 0; i < N; i++) {
            pos->y[i] = wrap_y(pos->y[i], box_ymin, Ly);
        }

        // Restore fixed atoms
        if (fixed_mask && fixed_copy) {
            for (int i = 0; i < N; i++) {
                if (fixed_mask[i]) {
                    pos->x[i] = fixed_copy->x[i];
                    pos->y[i] = fixed_copy->y[i];
                }
            }
        }
    }

    free(fx);
    free(fy);
    pa_free(fixed_copy);
}

/* ---------- cell-list based overlap deletion ---------- */
void iterative_delete_overlaps(PointArray *left, PointArray *right,
                               double threshold, double coord_cut,
                               double box_ymin, double box_ymax,
                               bool verbose) {
    double Ly = box_ymax - box_ymin;
    double xmin = 1e300, xmax = -1e300;

    // Compute x range from both arrays
    for (int i = 0; i < left->n; i++) {
        if (left->x[i] < xmin) xmin = left->x[i];
        if (left->x[i] > xmax) xmax = left->x[i];
    }
    for (int i = 0; i < right->n; i++) {
        if (right->x[i] < xmin) xmin = right->x[i];
        if (right->x[i] > xmax) xmax = right->x[i];
    }
    if (left->n == 0 || right->n == 0) return;
    double Lx = xmax - xmin;
    if (Lx < EPS) Lx = 1.0;

    // Cell size for overlap search
    double cell_size = threshold;
    if (cell_size < EPS) cell_size = 0.1;
    int nx = (int)ceil(Lx / cell_size);
    int ny = (int)ceil(Ly / cell_size);
    if (nx < 1) nx = 1;
    if (ny < 1) ny = 1;
    int ncell = nx * ny;

    bool changed;
    do {
        changed = false;
        int best_i = -1, best_j = -1;
        double best_dist = 1e9;

        // Allocate and build cell lists for left and right
        int *headL = malloc(ncell * sizeof(int));
        int *lsclL = malloc(left->n * sizeof(int));
        int *headR = malloc(ncell * sizeof(int));
        int *lsclR = malloc(right->n * sizeof(int));
        if (!headL || !lsclL || !headR || !lsclR) {
            fprintf(stderr, "Allocation error in cell-list deletion\n");
            exit(1);
        }
        for (int c = 0; c < ncell; c++) {
            headL[c] = -1;
            headR[c] = -1;
        }

        // Build left cell list
        for (int i = 0; i < left->n; i++) {
            int cx = (int)floor((left->x[i] - xmin) / cell_size);
            if (cx < 0) cx = 0;
            if (cx >= nx) cx = nx - 1;
            int cy = (int)floor((left->y[i] - box_ymin) / cell_size);
            cy = cy % ny;
            if (cy < 0) cy += ny;
            int cell = cy * nx + cx;
            lsclL[i] = headL[cell];
            headL[cell] = i;
        }

        // Build right cell list
        for (int i = 0; i < right->n; i++) {
            int cx = (int)floor((right->x[i] - xmin) / cell_size);
            if (cx < 0) cx = 0;
            if (cx >= nx) cx = nx - 1;
            int cy = (int)floor((right->y[i] - box_ymin) / cell_size);
            cy = cy % ny;
            if (cy < 0) cy += ny;
            int cell = cy * nx + cx;
            lsclR[i] = headR[cell];
            headR[cell] = i;
        }

        // Find closest overlapping pair using cell lists
        for (int i = 0; i < left->n; i++) {
            double xi = left->x[i];
            double yi = left->y[i];

            int cx = (int)floor((xi - xmin) / cell_size);
            if (cx < 0) cx = 0;
            if (cx >= nx) cx = nx - 1;
            int cy = (int)floor((yi - box_ymin) / cell_size);
            cy = cy % ny;
            if (cy < 0) cy += ny;

            for (int ii = cx - 1; ii <= cx + 1; ii++) {
                if (ii < 0 || ii >= nx) continue;  // no PBC in x
                for (int jj = cy - 1; jj <= cy + 1; jj++) {
                    int ccy = jj % ny;
                    if (ccy < 0) ccy += ny;
                    int cell = ccy * nx + ii;

                    for (int j = headR[cell]; j != -1; j = lsclR[j]) {
                        double dx = xi - right->x[j];
                        double dy = yi - right->y[j];
                        dy -= Ly * round(dy / Ly);
                        double dist = sqrt(dx*dx + dy*dy);
                        if (dist < threshold && dist < best_dist) {
                            best_dist = dist;
                            best_i = i;
                            best_j = j;
                        }
                    }
                }
            }
        }

        free(headL); free(lsclL);
        free(headR); free(lsclR);

        if (best_i < 0) break;  // no overlaps left

        // Determine which atom to delete based on coordination
        double lx = left->x[best_i], ly = left->y[best_i];
        double rx = right->x[best_j], ry = right->y[best_j];
        int coord_left = 0, coord_right = 0;

        for (int k = 0; k < left->n; k++) {
            if (k == best_i) continue;
            double dx = lx - left->x[k];
            double dy = ly - left->y[k];
            dy -= Ly * round(dy / Ly);
            if (sqrt(dx*dx + dy*dy) < coord_cut) coord_left++;
        }
        for (int k = 0; k < right->n; k++) {
            if (k == best_j) continue;
            double dx = rx - right->x[k];
            double dy = ry - right->y[k];
            dy -= Ly * round(dy / Ly);
            if (sqrt(dx*dx + dy*dy) < coord_cut) coord_right++;
        }

        if (coord_left >= coord_right) {
            // delete right atom
            for (int k = best_j; k < right->n - 1; k++) {
                right->x[k] = right->x[k+1];
                right->y[k] = right->y[k+1];
            }
            right->n--;
            if (verbose) printf("Deleted right atom, remaining right: %d\n", right->n);
        } else {
            // delete left atom
            for (int k = best_i; k < left->n - 1; k++) {
                left->x[k] = left->x[k+1];
                left->y[k] = left->y[k+1];
            }
            left->n--;
            if (verbose) printf("Deleted left atom, remaining left: %d\n", left->n);
        }
        changed = true;
    } while (changed);
}

/* ---------- process one translation ---------- */
typedef struct {
    PointArray *relaxed;
    double energy;
    PointArray *left_clean;
    PointArray *right_clean;
    int n_free;
} TranslationResult;

TranslationResult* process_translation(double dx, double dy,
                                       PointArray *left, PointArray *right,
                                       double threshold, double coord_cut,
                                       double box_xmin, double box_xmax,
                                       double box_ymin, double box_ymax,
                                       double rc, double sigma, double epsilon,
                                       double fixed_limit,
                                       bool verbose) {
    TranslationResult *res = malloc(sizeof(TranslationResult));
    if (!res) return NULL;

    double Ly = box_ymax - box_ymin;
    double half_dx = dx * 0.5;
    double half_dy = dy * 0.5;

    PointArray *shifted_left = pa_create(left->n);
    for (int i = 0; i < left->n; i++) {
        double x = left->x[i] - half_dx;
        double y = left->y[i] - half_dy;
        y = wrap_y(y, box_ymin, Ly);
        pa_append(shifted_left, x, y);
    }
    PointArray *shifted_right = pa_create(right->n);
    for (int i = 0; i < right->n; i++) {
        double x = right->x[i] + half_dx;
        double y = right->y[i] + half_dy;
        y = wrap_y(y, box_ymin, Ly);
        pa_append(shifted_right, x, y);
    }

    iterative_delete_overlaps(shifted_left, shifted_right, threshold, coord_cut,
                              box_ymin, box_ymax, verbose);

    PointArray *atoms = pa_create(shifted_left->n + shifted_right->n);
    for (int i = 0; i < shifted_left->n; i++)
        pa_append(atoms, shifted_left->x[i], shifted_left->y[i]);
    for (int i = 0; i < shifted_right->n; i++)
        pa_append(atoms, shifted_right->x[i], shifted_right->y[i]);

    PointArray *filtered = pa_create(atoms->n);
    for (int i = 0; i < atoms->n; i++) {
        if (fabs(atoms->x[i]) > EPS || fabs(atoms->y[i]) > EPS)
            pa_append(filtered, atoms->x[i], atoms->y[i]);
    }
    pa_free(atoms);
    atoms = filtered;

    if (atoms->n == 0) {
        res->relaxed = NULL;
        res->energy = 1e9;
        res->left_clean = NULL;
        res->right_clean = NULL;
        res->n_free = 0;
        pa_free(shifted_left); pa_free(shifted_right); pa_free(atoms);
        return res;
    }

    bool *fixed_mask = malloc(atoms->n * sizeof(bool));
    PointArray *fixed_pos = pa_create(atoms->n);
    int n_free = 0;
    for (int i = 0; i < atoms->n; i++) {
        if (fabs(atoms->x[i]) > fixed_limit) {
            fixed_mask[i] = true;
        } else {
            fixed_mask[i] = false;
            n_free++;
        }
        pa_append(fixed_pos, atoms->x[i], atoms->y[i]);
    }
    if (n_free == 0) {
        res->relaxed = NULL;
        res->energy = 1e9;
        res->left_clean = NULL;
        res->right_clean = NULL;
        res->n_free = 0;
        free(fixed_mask); pa_free(fixed_pos); pa_free(shifted_left);
        pa_free(shifted_right); pa_free(atoms);
        return res;
    }

    // Relaxation with faster, less strict settings
    gradient_descent_relax(atoms,
                           box_xmin, box_xmax, box_ymin, box_ymax,
                           rc, sigma, epsilon,
                           0.01, 0.0, 1500, 1e-3,
                           verbose,
                           fixed_mask, fixed_pos,
                           0.08 * 1.2);   // 8% of lattice constant max displacement

    double *fx = malloc(atoms->n * sizeof(double));
    double *fy = malloc(atoms->n * sizeof(double));
    double energy = lj_force_energy_cell(atoms, box_xmin, box_xmax,
                                         box_ymin, box_ymax, rc, sigma, epsilon,
                                         fx, fy);
    free(fx); free(fy);

    res->relaxed = atoms;
    res->energy = energy;
    res->left_clean = shifted_left;
    res->right_clean = shifted_right;
    res->n_free = n_free;

    free(fixed_mask);
    pa_free(fixed_pos);
    return res;
}

/* ---------- main ---------- */
int main() {
    double a = 1.2;   // lattice constant

    double sigma = a / pow(2.0, 1.0/6.0);
    double epsilon = 1.0;

    double box_xmin = -100.0 * a, box_xmax = 100.0 * a;
    double box_ymin = -36.373, box_ymax = 36.373;
    double margin = 3.0 * a;
    double rc = 3 * sigma;
    double r_eq = pow(2.0, 1.0/6.0) * sigma;
    double overlap_threshold = 0.5 * sigma;   // only delete extremely close atoms
    double coord_cut = 1.2 * r_eq;
    double fixed_limit = 95.0 * a;   // fix outermost ~5a layers

    double theta_deg = 25.0;
    double theta = theta_deg * PI / 180.0;

    // Scan up to 2*a with 0.2*a step
    double step = 0.2 * a;
    int n_dx = (int)(2.0 * a / step);
    int n_dy = n_dx;

    PointArray *lattice = generate_triangular_lattice(a, box_xmin, box_xmax,
                                                      box_ymin, box_ymax, margin);

    PointArray *left_all = pa_create(lattice->n);
    PointArray *right_all = pa_create(lattice->n);
    rotate_points(lattice, theta/2.0, left_all);
    rotate_points(lattice, -theta/2.0, right_all);

    PointArray *left_box = pa_create(left_all->n);
    PointArray *right_box = pa_create(right_all->n);
    filter_box(left_all, box_xmin, box_xmax, box_ymin, box_ymax, left_box);
    filter_box(right_all, box_xmin, box_xmax, box_ymin, box_ymax, right_box);

    PointArray *left = pa_create(left_box->n);
    PointArray *right = pa_create(right_box->n);
    for (int i = 0; i < left_box->n; i++) {
        if (left_box->x[i] < -EPS) pa_append(left, left_box->x[i], left_box->y[i]);
    }
    for (int i = 0; i < right_box->n; i++) {
        if (right_box->x[i] > EPS) pa_append(right, right_box->x[i], right_box->y[i]);
    }
    printf("Initial: left=%d, right=%d\n", left->n, right->n);

    // Best selection based on energy per free atom
    double best_per_free = 1e9;
    double best_energy = 1e9;
    double best_dx = 0, best_dy = 0;
    PointArray *best_relaxed = NULL;
    int best_free = 0;
    int total = (n_dx + 1) * (2*n_dy + 1);  // dx: -n_dx .. 0, dy: -n_dy .. +n_dy
    int count = 0;

    // Parallelize translations
    // dx scan restricted to non-positive values to avoid separating grains
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int ix = -n_dx; ix <= 0; ix++) {
        for (int iy = -n_dy; iy <= n_dy; iy++) {
            double dx = ix * step;
            double dy = iy * step;

            TranslationResult *res = process_translation(dx, dy,
                                                         left, right,
                                                         overlap_threshold, coord_cut,
                                                         box_xmin, box_xmax,
                                                         box_ymin, box_ymax,
                                                         rc, sigma, epsilon,
                                                         fixed_limit,
                                                         false);
            if (!res || !res->relaxed) {
                if (res) free(res);
                continue;
            }

            double energy = res->energy;
            int n_free = res->n_free;
            double per_free = energy / n_free;

            #pragma omp critical
            {
                count++;
                printf("Trying translation dx=%.3f, dy=%.3f (%d/%d)\n", dx, dy, count, total);
                printf("  Energy = %.6f (total), free atoms = %d, per free atom = %.6f\n",
                       energy, n_free, per_free);

                if (per_free < best_per_free) {
                    best_per_free = per_free;
                    best_energy = energy;
                    best_dx = dx; best_dy = dy;
                    best_free = n_free;

                    if (best_relaxed) pa_free(best_relaxed);
                    best_relaxed = pa_create(res->relaxed->n);
                    pa_copy(res->relaxed, best_relaxed);
                    printf("  *** New best per free atom energy: %.6f at dx=%.3f, dy=%.3f\n",
                           per_free, dx, dy);
                }
            }

            pa_free(res->relaxed);
            pa_free(res->left_clean);
            pa_free(res->right_clean);
            free(res);
        }
    }

    printf("\n=== Scan finished ===\n");
    printf("Best translation: dx=%.3f, dy=%.3f\n", best_dx, best_dy);
    if (best_relaxed) {
        printf("Best total energy: %.6f, free atoms: %d\n", best_energy, best_free);
        printf("Best energy per free atom: %.6f\n", best_per_free);

        PointArray *left_relaxed = pa_create(best_relaxed->n);
        PointArray *right_relaxed = pa_create(best_relaxed->n);
        for (int i = 0; i < best_relaxed->n; i++) {
            double x = best_relaxed->x[i], y = best_relaxed->y[i];
            if (x < -EPS) pa_append(left_relaxed, x, y);
            else if (x > EPS) pa_append(right_relaxed, x, y);
        }

        char fname[256];
        sprintf(fname, "bicrystal_scan_%.0fdeg_best.json", theta_deg);
        FILE *fp = fopen(fname, "w");
        if (fp) {
            fprintf(fp, "[\n");
            for (int i = 0; i < left_relaxed->n; i++) {
                fprintf(fp, "  {\"x\": %.12f, \"y\": %.12f, \"grain\": \"left\"}%s\n",
                        left_relaxed->x[i], left_relaxed->y[i],
                        (i == left_relaxed->n-1 && right_relaxed->n == 0) ? "" : ",");
            }
            for (int i = 0; i < right_relaxed->n; i++) {
                fprintf(fp, "  {\"x\": %.12f, \"y\": %.12f, \"grain\": \"right\"}%s\n",
                        right_relaxed->x[i], right_relaxed->y[i],
                        (i == right_relaxed->n-1) ? "" : ",");
            }
            fprintf(fp, "]\n");
            fclose(fp);
            printf("Best relaxed structure saved to %s\n", fname);
        }

        fp = fopen("best_structure.txt", "w");
        if (fp) {
            for (int i = 0; i < left_relaxed->n; i++)
                fprintf(fp, "%f %f 1\n", left_relaxed->x[i], left_relaxed->y[i]);
            for (int i = 0; i < right_relaxed->n; i++)
                fprintf(fp, "%f %f 2\n", right_relaxed->x[i], right_relaxed->y[i]);
            fclose(fp);
        }

        pa_free(left_relaxed);
        pa_free(right_relaxed);
        pa_free(best_relaxed);
    }

    pa_free(lattice);
    pa_free(left_all); pa_free(right_all);
    pa_free(left_box); pa_free(right_box);
    pa_free(left); pa_free(right);

    return 0;
}