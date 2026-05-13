#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>

/*
 * ObchodniCestujici8.c
 * ---------------------------------------------------------------------------
 * Travelling Salesman Problem (TSP) solved by Simulated Annealing (SA).
 *
 * Features:
 *  - Dimension switch: 2D or 3D (CLI parameter)
 *  - External temperature schedule file (step, temperature) with linear interpolation
 *  - If simulation step exceeds last schedule step, last temperature is kept
 *  - No T0/k parameters (temperature comes from the schedule)
 *  - For 2D mode, zsize is NOT required
 *
 * Output:
 *  - BASE.trip, BASE.set, BASE.conf, BASE.plb, BASE.fin, BASE.mol
 *  - .plb always contains 3 coordinates per city (z=0 in 2D), box[2]=0 in 2D
 *  - .conf/.fin are 2 columns in 2D; 3 columns in 3D
 *  - BASE.mol is overwritten if it exists
 *
 * Compile:
 *   gcc -O3 -march=native -ffast-math -funroll-loops -Wall -Wextra -pedantic -o ObchodniCestujici8_1 ObchodniCestujici8_1.c -lm
 */

#define POCET 6
#define MAX_BRUTE_FORCE_CITIES 12  /* brute-force allowed only for N <= this */

/* ------------------------------------------------------------------------- */
/* Helpers: I/O                                                              */
/* ------------------------------------------------------------------------- */

static int file_exists(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) return 0;
    fclose(f);
    return 1;
}

static FILE *xfopen(const char *path, const char *mode) {
    FILE *f = fopen(path, mode);
    if (!f) {
        perror(path);
        exit(EXIT_FAILURE);
    }
    return f;
}

/* ------------------------------------------------------------------------- */
/* Temperature schedule (linear interpolation)                               */
/* ------------------------------------------------------------------------- */

typedef struct {
    unsigned long step;
    float temp;
} temp_point_t;

typedef struct {
    unsigned long *steps;
    float *temps;
    size_t n;
    size_t idx; /* cached interval index for monotonic step queries */
} temp_schedule_t;

static int cmp_temp_point(const void *a, const void *b) {
    const temp_point_t *pa = (const temp_point_t *)a;
    const temp_point_t *pb = (const temp_point_t *)b;
    if (pa->step < pb->step) return -1;
    if (pa->step > pb->step) return 1;
    return 0;
}

static void free_temp_schedule(temp_schedule_t *s) {
    if (!s) return;
    free(s->steps);
    free(s->temps);
    s->steps = NULL;
    s->temps = NULL;
    s->n = 0;
    s->idx = 0;
}

static temp_schedule_t load_temp_schedule(const char *path) {
    if (!file_exists(path)) {
        fprintf(stderr, "Temperature schedule file does not exist: %s\n", path);
        exit(EXIT_FAILURE);
    }

    FILE *f = xfopen(path, "r");

    size_t cap = 128;
    size_t n = 0;
    temp_point_t *pts = (temp_point_t *)malloc(cap * sizeof(temp_point_t));
    if (!pts) {
        fprintf(stderr, "Out of memory while reading temperature schedule.\n");
        exit(EXIT_FAILURE);
    }

    /* Accept whitespace separated columns. Ignore trailing garbage lines. */
    while (1) {
        unsigned long step;
        float temp;
        int rc = fscanf(f, "%lu%f", &step, &temp);
        if (rc == EOF) break;
        if (rc != 2) {
            /* Skip the rest of the line */
            int c;
            while ((c = fgetc(f)) != '\n' && c != EOF) {}
            continue;
        }
        if (n == cap) {
            cap *= 2;
            temp_point_t *tmp = (temp_point_t *)realloc(pts, cap * sizeof(temp_point_t));
            if (!tmp) {
                free(pts);
                fprintf(stderr, "Out of memory while reading temperature schedule.\n");
                exit(EXIT_FAILURE);
            }
            pts = tmp;
        }
        pts[n].step = step;
        pts[n].temp = temp;
        n++;
    }
    fclose(f);

    if (n == 0) {
        free(pts);
        fprintf(stderr, "Temperature schedule file is empty or unreadable: %s\n", path);
        exit(EXIT_FAILURE);
    }

    /* Sort by step. */
    qsort(pts, n, sizeof(temp_point_t), cmp_temp_point);

    /* Deduplicate steps (keep last). */
    size_t m = 0;
    for (size_t i = 0; i < n; i++) {
        if (m == 0 || pts[i].step != pts[m - 1].step) {
            pts[m++] = pts[i];
        } else {
            pts[m - 1] = pts[i];
        }
    }

    temp_schedule_t s;
    s.steps = (unsigned long *)malloc(m * sizeof(unsigned long));
    s.temps = (float *)malloc(m * sizeof(float));
    s.n = m;
    s.idx = 0;

    if (!s.steps || !s.temps) {
        free(pts);
        free_temp_schedule(&s);
        fprintf(stderr, "Out of memory allocating temperature schedule arrays.\n");
        exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < m; i++) {
        s.steps[i] = pts[i].step;
        s.temps[i] = pts[i].temp;
    }

    free(pts);
    return s;
}

static inline float schedule_temp(temp_schedule_t *s, unsigned long step) {
    if (s->n == 1) return s->temps[0];

    if (step <= s->steps[0]) {
        s->idx = 0;
        return s->temps[0];
    }

    /* advance cached index */
    while (s->idx + 1 < s->n && step > s->steps[s->idx + 1]) {
        s->idx++;
    }

    if (s->idx + 1 >= s->n) {
        return s->temps[s->n - 1];
    }

    const unsigned long s0 = s->steps[s->idx];
    const unsigned long s1 = s->steps[s->idx + 1];
    const float t0 = s->temps[s->idx];
    const float t1 = s->temps[s->idx + 1];

    if (step >= s1) {
        /* exact boundary */
        return t1;
    }

    const double denom = (double)(s1 - s0);
    const double frac = (denom > 0.0) ? ((double)(step - s0) / denom) : 0.0;
    return (float)((1.0 - frac) * (double)t0 + frac * (double)t1);
}

/* ------------------------------------------------------------------------- */
/* Setup (.set)                                                              */
/* ------------------------------------------------------------------------- */

static void write_setup(const char *path,
                        int dim,
                        int n_cities,
                        unsigned long max_cycles,
                        int print_period,
                        int variant,
                        float xsize,
                        float ysize,
                        float zsize,
                        int mix) {
    FILE *f = xfopen(path, "wt");
    fprintf(f, "Dim= %d\n", dim);
    fprintf(f, "Cityes= %d\n", n_cities);
    fprintf(f, "Cyc= %lu\n", max_cycles);
    fprintf(f, "Pr= %d\n", print_period);
    fprintf(f, "var= %d\n", variant);
    fprintf(f, "x= %f\n", xsize);
    fprintf(f, "y= %f\n", ysize);
    if (dim == 3) fprintf(f, "z= %f\n", zsize);
    fprintf(f, "mix= %d\n", mix);
    fclose(f);
}

static void read_setup(const char *path,
                       int *dim,
                       int *n_cities,
                       unsigned long *max_cycles,
                       int *print_period,
                       int *variant,
                       float *xsize,
                       float *ysize,
                       float *zsize,
                       int *mix) {
    if (!file_exists(path)) {
        fprintf(stderr, "Soubor neexistuje: %s\n", path);
        exit(EXIT_FAILURE);
    }

    FILE *f = xfopen(path, "r");

    if (fscanf(f, "Dim= %d\n", dim) != 1) {
        fprintf(stderr, "Nelze přečíst Dim\n");
        exit(EXIT_FAILURE);
    }
    if (fscanf(f, "Cityes= %d\n", n_cities) != 1) {
        fprintf(stderr, "Nelze přečíst Cityes\n");
        exit(EXIT_FAILURE);
    }
    if (fscanf(f, "Cyc= %lu\n", max_cycles) != 1) {
        fprintf(stderr, "Nelze přečíst Cyc\n");
        exit(EXIT_FAILURE);
    }
    if (fscanf(f, "Pr= %d\n", print_period) != 1) {
        fprintf(stderr, "Nelze přečíst Pr\n");
        exit(EXIT_FAILURE);
    }
    if (fscanf(f, "var= %d\n", variant) != 1) {
        fprintf(stderr, "Nelze přečíst var\n");
        exit(EXIT_FAILURE);
    }
    if (fscanf(f, "x= %f\n", xsize) != 1) {
        fprintf(stderr, "Nelze přečíst x\n");
        exit(EXIT_FAILURE);
    }
    if (fscanf(f, "y= %f\n", ysize) != 1) {
        fprintf(stderr, "Nelze přečíst y\n");
        exit(EXIT_FAILURE);
    }

    *zsize = 0.0f;
    if (*dim == 3) {
        if (fscanf(f, "z= %f\n", zsize) != 1) {
            fprintf(stderr, "Nelze přečíst z\n");
            exit(EXIT_FAILURE);
        }
    }

    if (fscanf(f, "mix= %d\n", mix) != 1) {
        fprintf(stderr, "Nelze přečíst mix\n");
        exit(EXIT_FAILURE);
    }

    fclose(f);
}

/* ------------------------------------------------------------------------- */
/* Config (.conf)                                                            */
/* ------------------------------------------------------------------------- */

static void read_config2(const char *path, int n_last, float *x, float *y) {
    if (!file_exists(path)) {
        fprintf(stderr, "Soubor neexistuje: %s\n", path);
        exit(EXIT_FAILURE);
    }

    FILE *f = xfopen(path, "r");
    int i = 0;
    while (i <= n_last && fscanf(f, "%f%f", &x[i], &y[i]) == 2) {
        i++;
    }
    fclose(f);

    if (i != n_last + 1) {
        fprintf(stderr, "Varování: načteno %d bodů, očekáváno %d\n", i, n_last + 1);
    }
}

static void read_config3(const char *path, int n_last, float *x, float *y, float *z) {
    if (!file_exists(path)) {
        fprintf(stderr, "Soubor neexistuje: %s\n", path);
        exit(EXIT_FAILURE);
    }

    FILE *f = xfopen(path, "r");
    int i = 0;
    while (i <= n_last && fscanf(f, "%f%f%f", &x[i], &y[i], &z[i]) == 3) {
        i++;
    }
    fclose(f);

    if (i != n_last + 1) {
        fprintf(stderr, "Varování: načteno %d bodů, očekáváno %d\n", i, n_last + 1);
    }
}

static void write_config_from_tour2(const char *path,
                                   int n_last,
                                   const float *px,
                                   const float *py,
                                   const int *tour) {
    FILE *f = xfopen(path, "wt");
    for (int i = 0; i <= n_last; i++) {
        const int id = tour[i];
        fprintf(f, "%f\t%f\n", px[id], py[id]);
    }
    fclose(f);
}

static void write_config_from_tour3(const char *path,
                                   int n_last,
                                   const float *px,
                                   const float *py,
                                   const float *pz,
                                   const int *tour) {
    FILE *f = xfopen(path, "wt");
    for (int i = 0; i <= n_last; i++) {
        const int id = tour[i];
        fprintf(f, "%f\t%f\t%f\n", px[id], py[id], pz[id]);
    }
    fclose(f);
}

/* Overwrite BASE.mol */
static void write_mol_header_overwrite(const char *path,
                                       int n_last,
                                       unsigned long max_cycles,
                                       int print_period,
                                       int dim,
                                       const char *tempfile) {
    FILE *molf = xfopen(path, "w");
    const int n_atoms = n_last + 1;

    fprintf(molf, "! Soubor vygenerovan v ramci reseni problemu obchodniho cestujiciho.\n");
    fprintf(molf,
            "!Number of cities=%d, Steps=%lu, Printed every %d steps, Dim=%d, TempSchedule=%s\n",
            n_atoms, max_cycles, print_period, dim, tempfile ? tempfile : "(none)");
    fprintf(molf, "\nnumber_of_atoms = %d\n\n\n\n", n_atoms);
    fprintf(molf, "atoms\n");
    fprintf(molf, "! i  atom-id  a-type charge chir nbonds bound_atoms\n");

    fprintf(molf, "%3d %3d-H \t H      0    0      2     %d %d\n", 0, 0, 1, n_last);

    for (int i = 1; i < n_last; i++) {
        fprintf(molf, "%3d %3d-H \t H      0    0      2     %d %d\n", i, i, i - 1, i + 1);
    }

    fprintf(molf, "%3d %3d-H \t H      0    0      2     %d %d\n", n_last, n_last, n_last - 1, 0);

    fclose(molf);
}

static void write_plb_header(const char *path, const float header[2]) {
    FILE *plb = xfopen(path, "wb");
    fwrite(header, sizeof(float), 2, plb);
    fclose(plb);
}

/* ------------------------------------------------------------------------- */
/* PLB helpers (buffered frame writing)                                      */
/* ------------------------------------------------------------------------- */

/*
 * Write a single PLB frame in one fwrite() call (much faster than
 * writing x/y/z with three fwrite() calls per city).
 *
 * Buffer must have size at least (3 + 3*N) floats.
 *
 * - If pz is NULL, z=0 is written for all cities (2D mode).
 */
static inline void plb_write_frame_tour(FILE *plb,
                                       const float box[3],
                                       int N,
                                       const int *tour,
                                       const float *px,
                                       const float *py,
                                       const float *pz,
                                       float *buf) {
    buf[0] = box[0];
    buf[1] = box[1];
    buf[2] = box[2];
    float *out = buf + 3;
    for (int t = 0; t < N; t++) {
        const int id = tour[t];
        out[3*t + 0] = px[id];
        out[3*t + 1] = py[id];
        out[3*t + 2] = pz ? pz[id] : 0.0f;
    }
    (void)fwrite(buf, sizeof(float), (size_t)(3 + 3*N), plb);
}

/*
 * Write a PLB frame when coordinates are already in tour order arrays
 * X/Y/Z of length N. If Z is NULL, z=0 is written.
 */
static inline void plb_write_frame_order(FILE *plb,
                                        const float box[3],
                                        int N,
                                        const float *X,
                                        const float *Y,
                                        const float *Z,
                                        float *buf) {
    buf[0] = box[0];
    buf[1] = box[1];
    buf[2] = box[2];
    float *out = buf + 3;
    for (int i = 0; i < N; i++) {
        out[3*i + 0] = X[i];
        out[3*i + 1] = Y[i];
        out[3*i + 2] = Z ? Z[i] : 0.0f;
    }
    (void)fwrite(buf, sizeof(float), (size_t)(3 + 3*N), plb);
}


static void append_runtime(const char *trip_path, const char *label, double seconds) {
    FILE *f = xfopen(trip_path, "a");

    int days = 0, hours = 0, minutes = 0;
    double sec = seconds;

    if (sec >= 86400.0) {
        days = (int)(sec / 86400.0);
        sec -= days * 86400.0;
    }
    if (sec >= 3600.0) {
        hours = (int)(sec / 3600.0);
        sec -= hours * 3600.0;
    }
    if (sec >= 60.0) {
        minutes = (int)(sec / 60.0);
        sec -= minutes * 60.0;
    }

    if (days > 0) {
        fprintf(f, "!%s run time was %d days %d hours %d minutes %f sec\n", label, days, hours,
                minutes, sec);
        printf("%s run time was %d days %d hours %d minutes %f sec\n", label, days, hours, minutes,
               sec);
    } else if (hours > 0) {
        fprintf(f, "!%s run time was %d hours %d minutes %f sec\n", label, hours, minutes, sec);
        printf("%s run time was %d hours %d minutes %f sec\n", label, hours, minutes, sec);
    } else if (minutes > 0) {
        fprintf(f, "!%s run time was %d minutes %f sec\n", label, minutes, sec);
        printf("%s run time was %d minutes %f sec\n", label, minutes, sec);
    } else {
        fprintf(f, "!%s run time was %f sec\n", label, sec);
        printf("%s run time was %f sec\n", label, sec);
    }

    fclose(f);
}

/* ------------------------------------------------------------------------- */
/* Fast RNG                                                                  */
/* ------------------------------------------------------------------------- */

typedef struct {
    unsigned int s;
} rng32_t;

static inline unsigned int xorshift32(rng32_t *r) {
    unsigned int x = r->s;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    r->s = x;
    return x;
}

static inline float rng01(rng32_t *r) {
    return (float)(xorshift32(r) >> 8) * (1.0f / 16777216.0f);
}

static inline int rng_index(rng32_t *r, int n_last) {
    return (int)(xorshift32(r) % (unsigned int)(n_last + 1));
}

static inline float rng_range(rng32_t *r, float a, float b) {
    return rng01(r) * (b - a) + a;
}

static void shuffle_int(rng32_t *r, int *arr, int n) {
    for (int i = n - 1; i > 0; i--) {
        int j = (int)(xorshift32(r) % (unsigned int)(i + 1));
        int t = arr[i];
        arr[i] = arr[j];
        arr[j] = t;
    }
}

/* ------------------------------------------------------------------------- */
/* Distances + tour length                                                   */
/* ------------------------------------------------------------------------- */

static float *build_distance_matrix2(int N, const float *px, const float *py) {
    const size_t NN = (size_t)N * (size_t)N;
    if (NN > (size_t)20000 * (size_t)20000) {
        fprintf(stderr, "N=%d too large for full distance matrix.\n", N);
        exit(EXIT_FAILURE);
    }

    float *dist = (float *)malloc(NN * sizeof(float));
    if (!dist) {
        fprintf(stderr, "Nedostatek paměti pro matici vzdáleností (%zu floatů).\n", NN);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < N; i++) {
        dist[(size_t)i * (size_t)N + (size_t)i] = 0.0f;
        for (int j = i + 1; j < N; j++) {
            const float dx = px[i] - px[j];
            const float dy = py[i] - py[j];
            const float d = sqrtf(dx * dx + dy * dy);
            dist[(size_t)i * (size_t)N + (size_t)j] = d;
            dist[(size_t)j * (size_t)N + (size_t)i] = d;
        }
    }
    return dist;
}

static float *build_distance_matrix3(int N, const float *px, const float *py, const float *pz) {
    const size_t NN = (size_t)N * (size_t)N;
    if (NN > (size_t)20000 * (size_t)20000) {
        fprintf(stderr, "N=%d too large for full distance matrix.\n", N);
        exit(EXIT_FAILURE);
    }

    float *dist = (float *)malloc(NN * sizeof(float));
    if (!dist) {
        fprintf(stderr, "Nedostatek paměti pro matici vzdáleností (%zu floatů).\n", NN);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < N; i++) {
        dist[(size_t)i * (size_t)N + (size_t)i] = 0.0f;
        for (int j = i + 1; j < N; j++) {
            const float dx = px[i] - px[j];
            const float dy = py[i] - py[j];
            const float dz = pz[i] - pz[j];
            const float d = sqrtf(dx * dx + dy * dy + dz * dz);
            dist[(size_t)i * (size_t)N + (size_t)j] = d;
            dist[(size_t)j * (size_t)N + (size_t)i] = d;
        }
    }
    return dist;
}

static inline float dmat(const float *dist, int N, int a, int b) {
    return dist[(size_t)a * (size_t)N + (size_t)b];
}

static float tour_length_from_tour(const float *dist, int N, const int *tour) {
    float sum = 0.0f;
    for (int i = 0; i < N - 1; i++) {
        sum += dmat(dist, N, tour[i], tour[i + 1]);
    }
    sum += dmat(dist, N, tour[N - 1], tour[0]);
    return sum;
}

/* 2-opt delta in O(1) on a tour permutation. i<k, not adjacent and not wrap-adjacent. */
static inline float delta_2opt_tour(const float *dist, int N, const int *tour, int i, int k) {
    const int im1 = (i == 0) ? (N - 1) : (i - 1);
    const int kp1 = (k == N - 1) ? 0 : (k + 1);

    const int a = tour[im1];
    const int b = tour[i];
    const int c = tour[k];
    const int d = tour[kp1];

    const float old_len = dmat(dist, N, a, b) + dmat(dist, N, c, d);
    const float new_len = dmat(dist, N, a, c) + dmat(dist, N, b, d);
    return new_len - old_len;
}

static inline void apply_2opt_tour(int *tour, int i, int k) {
    while (i < k) {
        int t = tour[i];
        tour[i] = tour[k];
        tour[k] = t;
        i++;
        k--;
    }
}


/* ------------------------------------------------------------------------- */
/* Brute force (exact TSP) for small N                                       */
/* ------------------------------------------------------------------------- */

/* Next permutation (lexicographic) for int array a[0..n-1]. Returns 1 if permuted, 0 if wrapped. */
static int next_permutation_int(int *a, int n) {
    int i = n - 2;
    while (i >= 0 && a[i] >= a[i + 1]) i--;
    if (i < 0) return 0;
    int j = n - 1;
    while (a[j] <= a[i]) j--;
    int tmp = a[i]; a[i] = a[j]; a[j] = tmp;
    for (int l = i + 1, r = n - 1; l < r; l++, r--) {
        tmp = a[l]; a[l] = a[r]; a[r] = tmp;
    }
    return 1;
}

/* Exact TSP by enumerating (N-1)!/2 tours:
   - Fix city 0 as start
   - Enumerate permutations of cities 1..N-1
   - Skip reverse duplicates by enforcing perm[0] < perm[N-2]
*/
static float brute_force_best_tour(const float *dist, int N, int *best_tour_out) {
    const int M = N - 1;
    int *perm = (int *)malloc((size_t)M * sizeof(int));
    if (!perm) {
        fprintf(stderr, "Out of memory in brute force.\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < M; i++) perm[i] = i + 1; /* 1..N-1 */

    float best_len = INFINITY;
    int first = 1;
    do {
        /* remove reverse duplicates */
        if (perm[0] > perm[M - 1]) continue;

        float sum = 0.0f;
        int prev = 0;
        for (int i = 0; i < M; i++) {
            int cur = perm[i];
            sum += dmat(dist, N, prev, cur);
            prev = cur;
        }
        sum += dmat(dist, N, prev, 0);

        if (first || sum < best_len) {
            first = 0;
            best_len = sum;
            if (best_tour_out) {
                best_tour_out[0] = 0;
                for (int i = 0; i < M; i++) best_tour_out[i + 1] = perm[i];
            }
        }
    } while (next_permutation_int(perm, M));

    free(perm);
    return best_len;
}

/* ------------------------------------------------------------------------- */
/* "My correction" 2D/3D (kept, operating on coordinate order)              */
/* ------------------------------------------------------------------------- */

static int factorial(int n) {
    if (n <= 1) return 1;
    return n * factorial(n - 1);
}

static int wrap_index(int i, int n_last) {
    if (i >= 0) {
        return (i <= n_last) ? i : (i - n_last - 1);
    }
    return n_last + i + 1;
}

static void window_cities(int n_last, int city, int *a) {
    for (int i = 0; i < POCET; i++) {
        int v = city + i;
        if (v > n_last) v -= (n_last + 1);
        if (v < 0) v += (n_last + 1);
        a[i] = v;
    }
}

static int find_in_window(const int *arr, int value) {
    for (int i = 0; i < POCET; i++) {
        if (arr[i] == value) return i;
    }
    return -1;
}

static void nth_permutation_window(int n_last, int nth, int *D, int city) {
    const int fact = factorial(POCET);
    if (nth < 0 || nth >= fact) {
        nth %= fact;
        if (nth < 0) nth += fact;
    }

    int base[POCET];
    window_cities(n_last, city, base);

    int a[POCET];
    for (int i = 0; i < POCET; i++) a[i] = i;

    if (nth == 0) {
        for (int i = 0; i < POCET; i++) D[i] = base[i];
        return;
    }

    for (int step = 0; step < nth; step++) {
        int i = POCET - 2;
        while (i >= 0 && a[i] > a[i + 1]) i--;
        int j = POCET - 1;
        while (a[j] < a[i]) j--;
        int tmp = a[j];
        a[j] = a[i];
        a[i] = tmp;

        i++;
        for (j = POCET - 1; j > i; i++, j--) {
            tmp = a[i];
            a[i] = a[j];
            a[j] = tmp;
        }
    }

    for (int i = 0; i < POCET; i++) D[i] = base[a[i]];
}

/* ---- 2D correction ---- */

static float dist_xy_order(const float *x, const float *y, int i, int j) {
    const float dx = x[i] - x[j];
    const float dy = y[i] - y[j];
    return sqrtf(dx * dx + dy * dy);
}

static float tour_length_order2(const float *x, const float *y, int n_last) {
    float sum = 0.0f;
    for (int i = 0; i <= n_last; i++) {
        const int j = (i == n_last) ? 0 : (i + 1);
        sum += dist_xy_order(x, y, i, j);
    }
    return sum;
}

static float window_path_length2(const int *order,
                                 const float *x,
                                 const float *y,
                                 int city,
                                 int n_last) {
    int list[POCET + 2];
    list[0] = wrap_index(city - 1, n_last);
    list[POCET + 1] = wrap_index(city + POCET, n_last);

    for (int i = 0; i < POCET; i++) list[i + 1] = order[i];

    float sum = 0.0f;
    for (int i = 0; i < POCET + 1; i++) {
        sum += dist_xy_order(x, y, list[i], list[i + 1]);
    }
    return sum;
}

static void improve_city_window2(float *x, float *y, int city, int n_last, int variant) {
    const int fact = factorial(POCET);

    int best_order[POCET];
    float best_value = 0.0f;

    for (int i = 0; i < fact; i++) {
        int order[POCET];
        nth_permutation_window(n_last, i, order, city);

        float val = (float)variant * window_path_length2(order, x, y, city, n_last);
        if (i == 0 || val <= best_value) {
            best_value = val;
            memcpy(best_order, order, sizeof(best_order));
        }
    }

    int window_ids[POCET];
    for (int i = 0; i < POCET; i++) window_ids[i] = wrap_index(city + i, n_last);

    float x_tmp[POCET];
    float y_tmp[POCET];
    for (int i = 0; i < POCET; i++) {
        x_tmp[i] = x[window_ids[i]];
        y_tmp[i] = y[window_ids[i]];
    }

    for (int i = 0; i < POCET; i++) {
        const int idx_in_old = find_in_window(window_ids, best_order[i]);
        if (idx_in_old < 0) continue;
        x[window_ids[i]] = x_tmp[idx_in_old];
        y[window_ids[i]] = y_tmp[idx_in_old];
    }
}

static double my_correction2(float *x, float *y, int n_last, int variant) {
    clock_t start = clock();
    for (int i = 0; i <= n_last; i++) {
        improve_city_window2(x, y, i, n_last, variant);
    }
    clock_t finish = clock();
    return (double)(finish - start) / (double)CLOCKS_PER_SEC;
}

/* ---- 3D correction ---- */

static float dist_xyz_order(const float *x, const float *y, const float *z, int i, int j) {
    const float dx = x[i] - x[j];
    const float dy = y[i] - y[j];
    const float dz = z[i] - z[j];
    return sqrtf(dx * dx + dy * dy + dz * dz);
}

static float tour_length_order3(const float *x, const float *y, const float *z, int n_last) {
    float sum = 0.0f;
    for (int i = 0; i <= n_last; i++) {
        const int j = (i == n_last) ? 0 : (i + 1);
        sum += dist_xyz_order(x, y, z, i, j);
    }
    return sum;
}

static float window_path_length3(const int *order,
                                 const float *x,
                                 const float *y,
                                 const float *z,
                                 int city,
                                 int n_last) {
    int list[POCET + 2];
    list[0] = wrap_index(city - 1, n_last);
    list[POCET + 1] = wrap_index(city + POCET, n_last);

    for (int i = 0; i < POCET; i++) list[i + 1] = order[i];

    float sum = 0.0f;
    for (int i = 0; i < POCET + 1; i++) {
        sum += dist_xyz_order(x, y, z, list[i], list[i + 1]);
    }
    return sum;
}

static void improve_city_window3(float *x, float *y, float *z,
                                 int city, int n_last, int variant) {
    const int fact = factorial(POCET);

    int best_order[POCET];
    float best_value = 0.0f;

    for (int i = 0; i < fact; i++) {
        int order[POCET];
        nth_permutation_window(n_last, i, order, city);

        float val = (float)variant * window_path_length3(order, x, y, z, city, n_last);
        if (i == 0 || val <= best_value) {
            best_value = val;
            memcpy(best_order, order, sizeof(best_order));
        }
    }

    int window_ids[POCET];
    for (int i = 0; i < POCET; i++) window_ids[i] = wrap_index(city + i, n_last);

    float x_tmp[POCET], y_tmp[POCET], z_tmp[POCET];
    for (int i = 0; i < POCET; i++) {
        x_tmp[i] = x[window_ids[i]];
        y_tmp[i] = y[window_ids[i]];
        z_tmp[i] = z[window_ids[i]];
    }

    for (int i = 0; i < POCET; i++) {
        const int idx_in_old = find_in_window(window_ids, best_order[i]);
        if (idx_in_old < 0) continue;
        x[window_ids[i]] = x_tmp[idx_in_old];
        y[window_ids[i]] = y_tmp[idx_in_old];
        z[window_ids[i]] = z_tmp[idx_in_old];
    }
}

static double my_correction3(float *x, float *y, float *z, int n_last, int variant) {
    clock_t start = clock();
    for (int i = 0; i <= n_last; i++) {
        improve_city_window3(x, y, z, i, n_last, variant);
    }
    clock_t finish = clock();
    return (double)(finish - start) / (double)CLOCKS_PER_SEC;
}

/* ------------------------------------------------------------------------- */
/* Simulated Annealing (optimized)                                           */
/* ------------------------------------------------------------------------- */

static double run_simulated_annealing_fast(int *tour,
                                           int N,
                                           const float *px,
                                           const float *py,
                                           const float *pz,
                                           const float *dist,
                                           int variant,
                                           unsigned long max_cycles,
                                           int print_period,
                                           float xsize,
                                           float ysize,
                                           float zsize,
                                           const float box_size[3],
                                           const char *path_plb,
                                           const char *path_trip,
                                           float *temperature,
                                           unsigned long *cycle,
                                           float *current_len,
                                           rng32_t *rng,
                                           temp_schedule_t *schedule) {
    /*
     * Simulated annealing on a tour permutation with O(1) 2-opt delta.
     *
     * Performance notes:
     *  - temperature is obtained via cached linear interpolation (schedule->idx)
     *  - tour length is updated incrementally (no full recompute in the loop)
     *  - PLB output is buffered and written in one fwrite() call per frame
     */

    /* scale factor for acceptance; use box diagonal (zsize=0 in 2D) */
    const float korr = sqrtf(xsize * xsize + ysize * ysize + zsize * zsize);
    const float inv_korr = (korr > 0.0f) ? (1.0f / korr) : 1.0f;

    FILE *plb = xfopen(path_plb, "ab");
    FILE *trip = xfopen(path_trip, "a");

    /* Large buffered I/O reduces syscall overhead for long runs. */
    setvbuf(plb, NULL, _IOFBF, 1 << 20);
    setvbuf(trip, NULL, _IOFBF, 1 << 20);

    float *frame_buf = (float *)malloc((size_t)(3 + 3 * N) * sizeof(float));
    if (!frame_buf) {
        fprintf(stderr, "Out of memory allocating PLB frame buffer.\n");
        exit(EXIT_FAILURE);
    }

    clock_t start = clock();

    while (*cycle <= max_cycles) {
        *temperature = schedule_temp(schedule, *cycle);
        if (*temperature < 1.0e-12f) *temperature = 1.0e-12f;
        const float inv_Tk = inv_korr / (*temperature);

        /* Pick a valid 2-opt move (i,k): not adjacent and not wrap-adjacent. */
        int i, kpos;
        do {
            i = rng_index(rng, N - 1);
            kpos = rng_index(rng, N - 1);
            if (kpos < i) {
                int tmp = i;
                i = kpos;
                kpos = tmp;
            }
        } while ((kpos - i) < 2 || (i == 0 && kpos == N - 1));

        const float delta = delta_2opt_tour(dist, N, tour, i, kpos);
        const float dE = (float)variant * delta;

        int accept;
        if (dE < 0.0f) {
            accept = 1;
        } else {
            const float p = expf(-dE * inv_Tk);
            accept = (p > rng01(rng));
        }

        if (accept) {
            apply_2opt_tour(tour, i, kpos);
            *current_len += delta;
        }

        if ((*cycle) % (unsigned long)print_period == 0) {
            plb_write_frame_tour(plb, box_size, N, tour, px, py, pz, frame_buf);
            fprintf(trip, "%15lu\t%lu\t%f\t%f\n",
                    *cycle,
                    (unsigned long)(*cycle / (unsigned long)print_period),
                    *temperature,
                    fabsf((float)variant * (*current_len)));
        }

        (*cycle)++;
    }

    /* Ensure final state is printed even if max_cycles not aligned with print_period. */
    if (max_cycles % (unsigned long)print_period != 0) {
        *temperature = schedule_temp(schedule, max_cycles);
        if (*temperature < 1.0e-12f) *temperature = 1.0e-12f;

        plb_write_frame_tour(plb, box_size, N, tour, px, py, pz, frame_buf);
        fprintf(trip, "%15lu\t%lu\t%f\t%f\n",
                max_cycles,
                (unsigned long)(max_cycles / (unsigned long)print_period),
                *temperature,
                fabsf((float)variant * (*current_len)));
    }

    fflush(plb);
    fflush(trip);
    fclose(plb);
    fclose(trip);

    free(frame_buf);

    clock_t finish = clock();
    return (double)(finish - start) / (double)CLOCKS_PER_SEC;
}

/* ------------------------------------------------------------------------- */
/* CLI                                                                       */
/* ------------------------------------------------------------------------- */

static void print_usage(const char *prog) {
    printf("WRONG number of arguments.\n");
    printf("Program for solving Travelling salesman problem by simulated annealing (2D/3D).\n");
    printf("Number of cities must be higher than 4.\n\n");

    printf("Usage (generate random cities):\n");
    printf("  %s DIM N NAME Cyc Pr xsize ysize [zsize] TEMPFILE [maximize] [-f]\n", prog);
    printf("    DIM=2 or DIM=3\n");
    printf("    TEMPFILE: two columns: step temperature (linear interpolation)\n");
    printf("\n");

    printf("Usage (load from NAME.conf + NAME.set):\n");
    printf("  %s DIM NAME TEMPFILE [maximize] [-f]\n", prog);
    printf("\n");

    printf("Notes:\n");
    printf("  - If simulation is longer than the last step in TEMPFILE, the last temperature is kept.\n");
    printf("  - In 2D mode, zsize is not required; z is set to 0.\n");
    printf("  - BASE.mol, BASE.plb, BASE.trip will be overwritten if they exist!\n\n");

    printf("Flags:\n");
    printf("  maximize : maximize distance (default is minimize)\n");
    printf("  -f       : after SA+correction run brute-force exact search (only for small N)\n\n");

    printf("Example 2D:\n");
    printf("  %s 2 50 Vystup 5000000 20000 1 1 temp.dat\n", prog);
    printf("Example 3D:\n");
    printf("  %s 3 50 Vystup 5000000 20000 1 1 1 temp.dat\n", prog);
}

int main(int argc, char *argv[]) {
    int dim = 0;
    int N = 0;
    int print_period = 100;
    int variant = 1;   /* 1=minimize, -1=maximize */
    int mix = 1;

    int do_bruteforce = 0;

    unsigned long max_cycles = 5000;
    unsigned long cycle = 1;

    float xsize = 1.0f, ysize = 1.0f, zsize = 0.0f;
    float temperature = 0.0f;

    /* seed fast RNG */
    rng32_t rng;
    rng.s = (unsigned int)time(NULL) ^ (unsigned int)(uintptr_t)&rng;
    if (rng.s == 0) rng.s = 0xA341316Cu;

    /*
      Two modes:
        - generate:
            DIM N NAME Cyc Pr xsize ysize [zsize] TEMPFILE [maximize] [-f]
        - load:
            DIM NAME TEMPFILE [maximize] [-f]
    */

    if (argc < 4) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    dim = atoi(argv[1]);
    if (!(dim == 2 || dim == 3)) {
        fprintf(stderr, "DIM must be 2 or 3 (got %d)\n", dim);
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    int is_load_mode = 0;
    if (argc >= 4) {
        /* load mode: DIM NAME TEMPFILE ... (NAME not numeric) */
        /* We decide by argument count and the position of mandatory args. */
        /* If argv[2] is not an integer or argc matches load mode minimal counts, treat as load. */
        char *endp = NULL;
        long maybeN = strtol(argv[2], &endp, 10);
        (void)maybeN;
        if (endp == argv[2] || *endp != '\0') {
            is_load_mode = 1;
        } else {
            is_load_mode = 0;
        }
        /* However, "generate" always needs at least: prog DIM N NAME Cyc Pr x y TEMPFILE => 9 args (2D)
           or plus zsize => 10 args (3D). If too few, it's definitely load mode. */
        if ((!is_load_mode) && ((dim == 2 && argc < 9) || (dim == 3 && argc < 10))) {
            is_load_mode = 1;
        }
    }

    const char *tempfile = NULL;
    char base[255] = {0};

    int flags_start = 0;

    if (is_load_mode) {
        /* DIM NAME TEMPFILE [flags...] */
        if (argc < 4) {
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }
        snprintf(base, sizeof(base), "%s", argv[2]);
        tempfile = argv[3];
        flags_start = 4;
    } else {
        /* generate mode */
        N = atoi(argv[2]);
        snprintf(base, sizeof(base), "%s", argv[3]);
        max_cycles = strtoul(argv[4], NULL, 10);
        print_period = atoi(argv[5]);
        xsize = (float)atof(argv[6]);
        ysize = (float)atof(argv[7]);

        if (dim == 3) {
            if (argc < 10) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            zsize = (float)atof(argv[8]);
            tempfile = argv[9];
            flags_start = 10;
        } else {
            zsize = 0.0f;
            tempfile = argv[8];
            flags_start = 9;
        }
    }

    /* Parse flags: "maximize" and "-f" can appear in any order after required args. */
    for (int i = flags_start; i < argc; i++) {
        if (strcmp(argv[i], "maximize") == 0) {
            variant = -1;
        } else if (strcmp(argv[i], "-f") == 0) {
            do_bruteforce = 1;
        } else {
            fprintf(stderr, "Unknown argument: %s\n", argv[i]);
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }
    }

    char path_trip[300], path_conf[300], path_set[300], path_plb[300], path_mol[300], path_fin[300];
    snprintf(path_trip, sizeof(path_trip), "%s.trip", base);
    snprintf(path_conf, sizeof(path_conf), "%s.conf", base);
    snprintf(path_set, sizeof(path_set), "%s.set", base);
    snprintf(path_plb, sizeof(path_plb), "%s.plb", base);
    snprintf(path_mol, sizeof(path_mol), "%s.mol", base);
    snprintf(path_fin, sizeof(path_fin), "%s.fin", base);

    /* Temperature schedule is always external */
    temp_schedule_t schedule = load_temp_schedule(tempfile);

    if (is_load_mode) {
        int dim_from_set = 0;
        int variant_from_set = 1;
        read_setup(path_set, &dim_from_set, &N, &max_cycles, &print_period, &variant_from_set,
                   &xsize, &ysize, &zsize, &mix);

        if (dim_from_set != dim) {
            fprintf(stderr,
                    "Warning: DIM on CLI (%d) differs from Dim in %s (%d). Using CLI DIM.\n",
                    dim, path_set, dim_from_set);
        }

        /* If no maximize flag was given, keep variant from .set */
        if (variant != -1) variant = variant_from_set;

        /* In 2D mode, force zsize=0 even if set says otherwise. */
        if (dim == 2) zsize = 0.0f;
    }

    if (N <= 4) {
        fprintf(stderr, "Number of cities must be higher than 4 (N=%d)\n", N);
        free_temp_schedule(&schedule);
        return EXIT_FAILURE;
    }

    const int n_last = N - 1;

    /* If brute force requested but too large N, warn and disable */
    if (do_bruteforce && N > MAX_BRUTE_FORCE_CITIES) {
        fprintf(stderr,
                "Warning: -f requested but N=%d is larger than MAX_BRUTE_FORCE_CITIES=%d. Ignoring -f.\n",
                N, MAX_BRUTE_FORCE_CITIES);
        do_bruteforce = 0;
    }

    /* PLB: keep 3 coords always, but zsize=0 in 2D */
    float box_size[3] = {xsize, ysize, (dim == 3) ? zsize : 0.0f};
    float header[2] = {(float)N, -3.0f};

    float *px = (float *)malloc((size_t)N * sizeof(float));
    float *py = (float *)malloc((size_t)N * sizeof(float));
    float *pz = (float *)malloc((size_t)N * sizeof(float)); /* allocated even for 2D */
    int *tour = (int *)malloc((size_t)N * sizeof(int));
    if (!px || !py || !pz || !tour) {
        fprintf(stderr, "Nedostatek paměti\n");
        free(px); free(py); free(pz); free(tour);
        free_temp_schedule(&schedule);
        return EXIT_FAILURE;
    }

    if (is_load_mode) {
        if (dim == 2) {
            read_config2(path_conf, n_last, px, py);
            for (int i = 0; i < N; i++) pz[i] = 0.0f;
        } else {
            read_config3(path_conf, n_last, px, py, pz);
        }
    } else {
        for (int i = 0; i < N; i++) {
            px[i] = rng_range(&rng, 0.0f, xsize);
            py[i] = rng_range(&rng, 0.0f, ysize);
            pz[i] = (dim == 3) ? rng_range(&rng, 0.0f, zsize) : 0.0f;
        }
    }

    for (int i = 0; i < N; i++) tour[i] = i;
    if (mix != 0) shuffle_int(&rng, tour, N);

    /* Distances */
    float *dist = NULL;
    if (dim == 2) dist = build_distance_matrix2(N, px, py);
    else dist = build_distance_matrix3(N, px, py, pz);

    /* Write setup (temperature schedule file is not stored in .set) */
    write_setup(path_set, dim, N, max_cycles, print_period, variant, xsize, ysize, zsize, mix);

    /* Output headers (overwrite) */
    write_plb_header(path_plb, header);
    write_mol_header_overwrite(path_mol, n_last, max_cycles, print_period, dim, tempfile);

    /* Init .trip */
    {
        FILE *tripf = xfopen(path_trip, "w");
        fprintf(tripf,
                "!Number of cities=%d, Steps=%lu, Printed every %d steps, Dim=%d, TempSchedule=%s\n",
                N, max_cycles, print_period, dim, tempfile);
        fprintf(tripf, "!Steps\tFrames\tTemperature\tLength of trip\n");
        fclose(tripf);
    }

    /* Initial length */
    float current_len = tour_length_from_tour(dist, N, tour);

    /* Initial temperature and initial frame */
    temperature = schedule_temp(&schedule, cycle);
    if (temperature < 1.0e-12f) temperature = 1.0e-12f;

    {
        FILE *plb = xfopen(path_plb, "ab");
        setvbuf(plb, NULL, _IOFBF, 1 << 20);
        float *frame_buf = (float *)malloc((size_t)(3 + 3 * N) * sizeof(float));
        if (!frame_buf) {
            fprintf(stderr, "Out of memory allocating PLB frame buffer.\n");
            exit(EXIT_FAILURE);
        }
        plb_write_frame_tour(plb, box_size, N, tour, px, py, pz, frame_buf);
        free(frame_buf);
        fclose(plb);
    }

    {
        FILE *tripf = xfopen(path_trip, "a");
        fprintf(tripf, "%15lu\t%lu\t%f\t%f\n",
                cycle,
                (unsigned long)(cycle / (unsigned long)print_period),
                temperature,
                fabsf((float)variant * current_len));
        fclose(tripf);
    }

    /* SA */
    const double sa_time = run_simulated_annealing_fast(tour,
                                                        N,
                                                        px,
                                                        py,
                                                        pz,
                                                        dist,
                                                        variant,
                                                        max_cycles,
                                                        print_period,
                                                        xsize,
                                                        ysize,
                                                        (dim == 3) ? zsize : 0.0f,
                                                        box_size,
                                                        path_plb,
                                                        path_trip,
                                                        &temperature,
                                                        &cycle,
                                                        &current_len,
                                                        &rng,
                                                        &schedule);

    append_runtime(path_trip, "Simulated annealing", sa_time);
    printf("Final trip length from simulated annealing is: %f\n", current_len);

    /* Materialize tour order for correction */
    float *X = (float *)malloc((size_t)N * sizeof(float));
    float *Y = (float *)malloc((size_t)N * sizeof(float));
    float *Z = (float *)malloc((size_t)N * sizeof(float));
    if (!X || !Y || !Z) {
        fprintf(stderr, "Nedostatek paměti\n");
        free(dist); free(px); free(py); free(pz); free(tour);
        free(X); free(Y); free(Z);
        free_temp_schedule(&schedule);
        return EXIT_FAILURE;
    }

    for (int i = 0; i < N; i++) {
        const int id = tour[i];
        X[i] = px[id];
        Y[i] = py[id];
        Z[i] = (dim == 3) ? pz[id] : 0.0f;
    }

    /* My correction */
    double corr_time_total = 0.0;
    float final_sum = current_len;

    if (N <= POCET) {
        fprintf(stderr, "WARNING: N=%d <= POCET=%d, skipping My correction.\n", N, POCET);
        FILE *tripf = xfopen(path_trip, "a");
        fprintf(tripf, "!My correction skipped because N=%d <= POCET=%d\n", N, POCET);
        fclose(tripf);
    } else {
        for (int r = 0; r <= 3; r++) {
        if (dim == 2) {
            corr_time_total += my_correction2(X, Y, n_last, variant);
            final_sum = tour_length_order2(X, Y, n_last);
        } else {
            corr_time_total += my_correction3(X, Y, Z, n_last, variant);
            final_sum = tour_length_order3(X, Y, Z, n_last);
        }

        FILE *plb = xfopen(path_plb, "ab");
        setvbuf(plb, NULL, _IOFBF, 1 << 20);
        float *frame_buf = (float *)malloc((size_t)(3 + 3 * N) * sizeof(float));
        if (!frame_buf) {
            fprintf(stderr, "Out of memory allocating PLB frame buffer.\n");
            exit(EXIT_FAILURE);
        }
        plb_write_frame_order(plb, box_size, N, X, Y, (dim == 3) ? Z : NULL, frame_buf);
        free(frame_buf);
        fclose(plb);

        FILE *tripf = xfopen(path_trip, "a");
        fprintf(tripf, "%15lu\t%lu\t%f\t%f\n",
                cycle,
                (unsigned long)(cycle / (unsigned long)print_period),
                temperature,
                fabsf(final_sum));
        fclose(tripf);
    }

    }

    if (N > POCET) {
        append_runtime(path_trip, "My correction", corr_time_total);
        {
            FILE *tripf = xfopen(path_trip, "a");
            fprintf(tripf, "!Trip length after my correction follows:\n");
            fclose(tripf);
        }
        printf("Final length after My correction is: %f\n", final_sum);
    } else {
        printf("Final length (My correction skipped) is: %f\n", final_sum);
    }

    /* Write outputs */
    if (dim == 2) {
        write_config_from_tour2(path_conf, n_last, px, py, tour);
        FILE *f = xfopen(path_fin, "wt");
        for (int i = 0; i <= n_last; i++) fprintf(f, "%f\t%f\n", X[i], Y[i]);
        fclose(f);
    } else {
        write_config_from_tour3(path_conf, n_last, px, py, pz, tour);
        FILE *f = xfopen(path_fin, "wt");
        for (int i = 0; i <= n_last; i++) fprintf(f, "%f\t%f\t%f\n", X[i], Y[i], Z[i]);
        fclose(f);
    }

    /* Optional brute force exact search (-f) */
    if (do_bruteforce) {
        int *best_tour = (int *)malloc((size_t)N * sizeof(int));
        if (!best_tour) {
            fprintf(stderr, "Out of memory for brute force tour.\n");
        } else {
            clock_t bf_start = clock();
            const float best_len = brute_force_best_tour(dist, N, best_tour);
            clock_t bf_finish = clock();
            const double bf_time = (double)(bf_finish - bf_start) / (double)CLOCKS_PER_SEC;
            append_runtime(path_trip, "Brute force", bf_time);

            /* Append last frame to .plb (best tour order) */
            FILE *plb = xfopen(path_plb, "ab");
            setvbuf(plb, NULL, _IOFBF, 1 << 20);
            float *frame_buf = (float *)malloc((size_t)(3 + 3 * N) * sizeof(float));
            if (!frame_buf) {
                fprintf(stderr, "Out of memory allocating PLB frame buffer.\n");
                exit(EXIT_FAILURE);
            }
            plb_write_frame_tour(plb, box_size, N, best_tour, px, py, (dim == 3) ? pz : NULL, frame_buf);
            free(frame_buf);
            fclose(plb);

            /* Append to .trip */
            FILE *tripf = xfopen(path_trip, "a");
            fprintf(tripf, "!Best trip length:\n");
            fprintf(tripf, "%f\n", best_len);
            fclose(tripf);

            printf("Best trip length (brute force) is: %f\n", best_len);
            free(best_tour);
        }
    }

    free(dist);
    free(px);
    free(py);
    free(pz);
    free(tour);
    free(X);
    free(Y);
    free(Z);
    free_temp_schedule(&schedule);

    return EXIT_SUCCESS;
}
