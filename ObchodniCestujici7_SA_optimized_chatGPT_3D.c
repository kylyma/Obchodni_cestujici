#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>

/*
 * ObchodniCestujici7_SA_optimized_chatGPT_3D.c
 * ---------------------------------------------------------------------------
 * Travelling Salesman Problem (TSP) solved by Simulated Annealing (SA).
 *
 * 3D version:
 *  - Cities are points in 3D: (x,y,z)
 *  - Distances are Euclidean in 3D
 *  - .conf and .fin contain 3 columns (x y z)
 *  - .plb stores per-city coordinates (x,y,z) per frame
 *  - "My correction" local-window optimization works in 3D too
 *
 * Change (requested):
 *  - Output mol file is now OUTPUT_NAME.mol (i.e., BASE.mol)
 *  - BASE.mol is overwritten if it exists
 *
 * Compile (recommended):
 *   gcc -O3 -march=native -ffast-math -funroll-loops -Wall -Wextra -pedantic \
 *       -o ObchodniCestujici7_SA_optimized_chatGPT_3D ObchodniCestujici7_SA_optimized_chatGPT_3D.c -lm
 */

#define POCET 6

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

static void write_setup(const char *path,
                        int n_cities,
                        unsigned long max_cycles,
                        int print_period,
                        float t0,
                        float k,
                        int variant,
                        float xsize,
                        float ysize,
                        float zsize,
                        int mix) {
    FILE *f = xfopen(path, "wt");
    fprintf(f, "Cityes= %d\n", n_cities);
    fprintf(f, "Cyc= %lu\n", max_cycles);
    fprintf(f, "Pr= %d\n", print_period);
    fprintf(f, "T0= %f\n", t0);
    fprintf(f, "k= %f\n", k);
    fprintf(f, "var= %d\n", variant);
    fprintf(f, "x= %f\n", xsize);
    fprintf(f, "y= %f\n", ysize);
    fprintf(f, "z= %f\n", zsize);
    fprintf(f, "mix= %d\n", mix);
    fclose(f);
}

static void read_setup(const char *path,
                       int *n_cities,
                       unsigned long *max_cycles,
                       int *print_period,
                       float *t0,
                       float *k,
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
    if (fscanf(f, "T0= %f\n", t0) != 1) {
        fprintf(stderr, "Nelze přečíst T0\n");
        exit(EXIT_FAILURE);
    }
    if (fscanf(f, "k= %f\n", k) != 1) {
        fprintf(stderr, "Nelze přečíst k\n");
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
    if (fscanf(f, "z= %f\n", zsize) != 1) {
        fprintf(stderr, "Nelze přečíst z\n");
        exit(EXIT_FAILURE);
    }
    if (fscanf(f, "mix= %d\n", mix) != 1) {
        fprintf(stderr, "Nelze přečíst mix\n");
        exit(EXIT_FAILURE);
    }

    fclose(f);
}

/* Reads coordinates from a .conf file: x y z per line */
static void read_config3(const char *path, int n_last, float *x, float *y, float *z) {
    if (!file_exists(path)) {
        fprintf(stderr, "Soubor neexistuje: %s\n", path);
        exit(EXIT_FAILURE);
    }

    FILE *f = xfopen(path, "r");
    int i = 0;
    while (i <= n_last && fscanf(f, "%f\t%f\t%f", &x[i], &y[i], &z[i]) == 3) {
        i++;
    }
    fclose(f);

    if (i != n_last + 1) {
        fprintf(stderr, "Varování: načteno %d bodů, očekáváno %d\n", i, n_last + 1);
    }
}

/* Writes coordinates in tour order to a text file: x y z */
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

/* NOTE: This function now OVERWRITES the mol file (mode "w"). */
static void write_mol_header_overwrite(const char *path,
                                      int n_last,
                                      unsigned long max_cycles,
                                      int print_period,
                                      float t0,
                                      float k) {
    FILE *molf = xfopen(path, "w"); /* overwrite */
    const int n_atoms = n_last + 1;

    fprintf(molf, "! Soubor vygenerovan v ramci reseni problemu obchodniho cestujiciho.\n");
    fprintf(molf,
            "!Number of cities=%d, Steps=%lu, Printed every %d steps, Initial Temperature= %f, k=%f\n",
            n_atoms,
            max_cycles,
            print_period,
            t0,
            k);
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
/* Distances + tour length (3D)                                              */
/* ------------------------------------------------------------------------- */

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
/* "My correction" (3D, operating on coordinate order)                       */
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
/* Simulated Annealing (optimized, 3D output)                                */
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
                                           float t0,
                                           float k,
                                           float xsize,
                                           float ysize,
                                           float zsize,
                                           const float box_size[3],
                                           const char *path_plb,
                                           const char *path_trip,
                                           float *temperature,
                                           unsigned long *cycle,
                                           float *current_len,
                                           rng32_t *rng) {
    const float korr = sqrtf(xsize * xsize + ysize * ysize + zsize * zsize);
    const float inv_korr = (korr > 0.0f) ? (1.0f / korr) : 1.0f;

    FILE *plb = xfopen(path_plb, "ab");
    FILE *trip = xfopen(path_trip, "a");

    clock_t start = clock();

    while (*cycle <= max_cycles) {
        *temperature = max_cycles * k * t0 / (*cycle + (max_cycles * k - 1.0f));
        const float inv_Tk = inv_korr / (*temperature);

        int i = rng_index(rng, N - 1);
        int kpos = i;
        while (abs(i - kpos) < 2 || abs(i - kpos) == (N - 1)) {
            kpos = rng_index(rng, N - 1);
        }
        if (kpos < i) {
            int tmp = i;
            i = kpos;
            kpos = tmp;
        }

        const float delta = delta_2opt_tour(dist, N, tour, i, kpos);
        const float dE = (float)variant * delta;

        int accept = 0;
        if (dE < 0.0f) {
            accept = 1;
        } else {
            const float p = expf(-dE * inv_Tk);
            if (p > rng01(rng)) accept = 1;
        }

        if (accept) {
            apply_2opt_tour(tour, i, kpos);
            *current_len += delta;
        }

        if ((*cycle) % (unsigned long)print_period == 0) {
            fwrite(box_size, sizeof(float), 3, plb);
            for (int t = 0; t < N; t++) {
                const int id = tour[t];
                const float x = px[id];
                const float y = py[id];
                const float z = pz[id];
                fwrite(&x, sizeof(float), 1, plb);
                fwrite(&y, sizeof(float), 1, plb);
                fwrite(&z, sizeof(float), 1, plb);
            }

            fprintf(trip, "%15lu\t%lu\t%f\t%f\n",
                    *cycle,
                    (unsigned long)(*cycle / (unsigned long)print_period),
                    *temperature,
                    fabsf((float)variant * (*current_len)));
        }

        (*cycle)++;
    }

    if (max_cycles % (unsigned long)print_period != 0) {
        fwrite(box_size, sizeof(float), 3, plb);
        for (int t = 0; t < N; t++) {
            const int id = tour[t];
            const float x = px[id];
            const float y = py[id];
            const float z = pz[id];
            fwrite(&x, sizeof(float), 1, plb);
            fwrite(&y, sizeof(float), 1, plb);
            fwrite(&z, sizeof(float), 1, plb);
        }

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

    clock_t finish = clock();
    return (double)(finish - start) / (double)CLOCKS_PER_SEC;
}

/* ------------------------------------------------------------------------- */
/* CLI                                                                       */
/* ------------------------------------------------------------------------- */

static void print_usage(const char *prog) {
    printf("WRONG number of arguments.\n");
    printf("Program for solving Travelling salesman problem by simulated annealing (3D).\n");
    printf("Number of cities must be higher than 4.\n\n");
    printf("Usage:\n");
    printf("  %s N NAME Cyc Pr T0 k xsize ysize zsize [maximize]\n", prog);
    printf("  If optional last argument is present, distance will be MAXIMIZED.\n\n");
    printf("  N      number of cities\n");
    printf("  NAME   base name of output files\n");
    printf("  Cyc    number of cycles\n");
    printf("  Pr     print period in cycles\n");
    printf("  T0     initial temperature\n");
    printf("  k      small constant (e.g. 0.0001), larger => slower cooling\n");
    printf("  xsize, ysize, zsize  size of simulation box\n\n");
    printf("Example:\n");
    printf("  %s 50 Vystup 5000000 20000 200 0.0001 1 1 1\n\n", prog);
    printf("Alternative:\n");
    printf("  %s NAME\n", prog);
    printf("  Files NAME.conf and NAME.set must exist.\n\n");
    printf("WARNING:\n");
    printf("  NAME.plb NAME.trip and NAME.mol will be overwritten if exists!\n\n");
}

int main(int argc, char *argv[]) {
    int N = 0;
    int print_period = 100;
    int variant = 1;
    int mix = 1;

    unsigned long max_cycles = 5000;
    unsigned long cycle = 1;

    float xsize = 1.0f, ysize = 1.0f, zsize = 1.0f;
    float t0 = 0.0f, k = 0.0f;
    float temperature = 0.0f;

    rng32_t rng;
    rng.s = (unsigned int)time(NULL) ^ (unsigned int)(uintptr_t)&rng;
    if (rng.s == 0) rng.s = 0xA341316Cu;

    if ((argc != 10) && (argc != 11) && (argc != 2)) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    char base[255];
    snprintf(base, sizeof(base), "%s", (argc == 2) ? argv[1] : argv[2]);

    char path_trip[300], path_conf[300], path_set[300], path_plb[300], path_mol[300], path_fin[300];
    snprintf(path_trip, sizeof(path_trip), "%s.trip", base);
    snprintf(path_conf, sizeof(path_conf), "%s.conf", base);
    snprintf(path_set, sizeof(path_set), "%s.set", base);
    snprintf(path_plb, sizeof(path_plb), "%s.plb", base);

    /* CHANGED: write to BASE.mol (and overwrite) */
    snprintf(path_mol, sizeof(path_mol), "%s.mol", base);

    snprintf(path_fin, sizeof(path_fin), "%s.fin", base);

    if (argc == 2) {
        read_setup(path_set, &N, &max_cycles, &print_period, &t0, &k, &variant,
                   &xsize, &ysize, &zsize, &mix);
    } else {
        sscanf(argv[1], "%d", &N);
        sscanf(argv[3], "%lu", &max_cycles);
        sscanf(argv[4], "%d", &print_period);
        sscanf(argv[5], "%f", &t0);
        sscanf(argv[6], "%f", &k);
        sscanf(argv[7], "%f", &xsize);
        sscanf(argv[8], "%f", &ysize);
        sscanf(argv[9], "%f", &zsize);
        if (argc == 11) variant = -1;
    }

    if (N <= 4) {
        fprintf(stderr, "Number of cities must be higher than 4 (N=%d)\n", N);
        return EXIT_FAILURE;
    }

    const int n_last = N - 1;
    temperature = t0;

    float box_size[3] = {xsize, ysize, zsize};
    float header[2] = {(float)N, -3.0f};

    float *px = (float *)malloc((size_t)N * sizeof(float));
    float *py = (float *)malloc((size_t)N * sizeof(float));
    float *pz = (float *)malloc((size_t)N * sizeof(float));
    int *tour = (int *)malloc((size_t)N * sizeof(int));
    if (!px || !py || !pz || !tour) {
        fprintf(stderr, "Nedostatek paměti\n");
        return EXIT_FAILURE;
    }

    if (argc == 2) {
        read_config3(path_conf, n_last, px, py, pz);
    } else {
        for (int i = 0; i < N; i++) {
            px[i] = rng_range(&rng, 0.0f, xsize);
            py[i] = rng_range(&rng, 0.0f, ysize);
            pz[i] = rng_range(&rng, 0.0f, zsize);
        }
    }

    for (int i = 0; i < N; i++) tour[i] = i;
    if (mix != 0) shuffle_int(&rng, tour, N);

    float *dist = build_distance_matrix3(N, px, py, pz);

    write_setup(path_set, N, max_cycles, print_period, t0, k, variant, xsize, ysize, zsize, mix);

    write_plb_header(path_plb, header);

    /* CHANGED: overwrite BASE.mol */
    write_mol_header_overwrite(path_mol, n_last, max_cycles, print_period, t0, k);

    {
        FILE *tripf = xfopen(path_trip, "w");
        fprintf(tripf,
                "!Number of cities=%d, Steps=%lu, Printed every %d steps, Initial Temperature= %f, k=%f\n",
                N,
                max_cycles,
                print_period,
                t0,
                k);
        fprintf(tripf, "!Steps\tFrames\tTemperature\tLength of trip\n");
        fclose(tripf);
    }

    float current_len = tour_length_from_tour(dist, N, tour);

    {
        FILE *plb = xfopen(path_plb, "ab");
        fwrite(box_size, sizeof(float), 3, plb);
        for (int t = 0; t < N; t++) {
            const int id = tour[t];
            const float x = px[id];
            const float y = py[id];
            const float z = pz[id];
            fwrite(&x, sizeof(float), 1, plb);
            fwrite(&y, sizeof(float), 1, plb);
            fwrite(&z, sizeof(float), 1, plb);
        }
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

    const double sa_time = run_simulated_annealing_fast(tour,
                                                        N,
                                                        px,
                                                        py,
                                                        pz,
                                                        dist,
                                                        variant,
                                                        max_cycles,
                                                        print_period,
                                                        t0,
                                                        k,
                                                        xsize,
                                                        ysize,
                                                        zsize,
                                                        box_size,
                                                        path_plb,
                                                        path_trip,
                                                        &temperature,
                                                        &cycle,
                                                        &current_len,
                                                        &rng);

    append_runtime(path_trip, "Simulated annealing", sa_time);
    printf("Final trip length from simulated annealing is: %f\n", current_len);

    float *X = (float *)malloc((size_t)N * sizeof(float));
    float *Y = (float *)malloc((size_t)N * sizeof(float));
    float *Z = (float *)malloc((size_t)N * sizeof(float));
    if (!X || !Y || !Z) {
        fprintf(stderr, "Nedostatek paměti\n");
        return EXIT_FAILURE;
    }
    for (int i = 0; i < N; i++) {
        const int id = tour[i];
        X[i] = px[id];
        Y[i] = py[id];
        Z[i] = pz[id];
    }

    double corr_time_total = 0.0;
    float final_sum = current_len;

    for (int r = 0; r <= 3; r++) {
        corr_time_total += my_correction3(X, Y, Z, n_last, variant);
        final_sum = tour_length_order3(X, Y, Z, n_last);

        FILE *plb = xfopen(path_plb, "ab");
        fwrite(box_size, sizeof(float), 3, plb);
        for (int i = 0; i <= n_last; i++) {
            fwrite(&X[i], sizeof(float), 1, plb);
            fwrite(&Y[i], sizeof(float), 1, plb);
            fwrite(&Z[i], sizeof(float), 1, plb);
        }
        fclose(plb);

        FILE *tripf = xfopen(path_trip, "a");
        fprintf(tripf, "%15lu\t%lu\t%f\t%f\n",
                cycle,
                (unsigned long)(cycle / (unsigned long)print_period),
                temperature,
                fabsf(final_sum));
        fclose(tripf);
    }

    append_runtime(path_trip, "My correction", corr_time_total);
    {
        FILE *tripf = xfopen(path_trip, "a");
        fprintf(tripf, "!Trip length after my correction follows:\n");
        fclose(tripf);
    }
    printf("Final length after My correction is: %f\n", final_sum);

    write_config_from_tour3(path_conf, n_last, px, py, pz, tour);

    {
        FILE *f = xfopen(path_fin, "wt");
        for (int i = 0; i <= n_last; i++) fprintf(f, "%f\t%f\t%f\n", X[i], Y[i], Z[i]);
        fclose(f);
    }

    free(dist);
    free(px);
    free(py);
    free(pz);
    free(tour);
    free(X);
    free(Y);
    free(Z);

    return EXIT_SUCCESS;
}
