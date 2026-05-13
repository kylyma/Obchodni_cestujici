#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>

/*
 * ObchodniCestujici7_SA_ultrafast_A.c
 * ------------------------------------------------------------
 * Řešení problému obchodního cestujícího (TSP) metodou simulovaného žíhání.
 *
 * Poznámky k formátu dat:
 *  - Konfigurace je uložena jako pořadí bodů (měst) v polích X[], Y[].
 *  - "Cesta" je cyklus: poslední město se vrací do prvního.
 *
 * Tento soubor je upravená/verzovaná varianta původního kódu:
 *  - sjednocené formátování
 *  - doplněné komentáře k funkcím
 *  - zpřehledněné názvy proměnných a práce se soubory
 *  - opravené časté chyby (sizeof, malloc, strcat na argumentech atd.)
 *
 * Kompilace:
 *   gcc -O3 -Wall -Wextra -pedantic -o ObchodniCestujici7_SA_ultrafast_A ObchodniCestujici7_SA_ultrafast_A.c -lm
 *   gcc -O3 -march=native -ffast-math -funroll-loops -Wall -Wextra -pedantic -o tsp_ultraA ObchodniCestujici7_SA_ultrafast_A.c -lm
 */

/* Počet měst lokálně optimalizovaných v "MojeOprava" (okno pro permutace). */
#define POCET 6

/* ------------------------------------------------------------------------- */
/* Pomocné I/O a validace                                                    */
/* ------------------------------------------------------------------------- */

/* Vrátí 1, pokud soubor existuje a lze ho otevřít pro čtení, jinak 0. */
static int file_exists(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) {
        return 0;
    }
    fclose(f);
    return 1;
}

/* Bezpečné otevření souboru – při chybě ukončí program s hláškou. */
static FILE *xfopen(const char *path, const char *mode) {
    FILE *f = fopen(path, mode);
    if (!f) {
        perror(path);
        exit(EXIT_FAILURE);
    }
    return f;
}

/* Zapíše konfiguraci měst (pořadí) do textového souboru. */
static void write_config(const char *path, int n_last, const float *x, const float *y) {
    FILE *f = xfopen(path, "wt");
    for (int i = 0; i <= n_last; i++) {
        fprintf(f, "%f\t%f\n", x[i], y[i]);
    }
    fclose(f);
}

/* Zapíše parametrický soubor .set. */
static void write_setup(const char *path,
                        int n_cities,
                        unsigned long max_cycles,
                        int print_period,
                        float t0,
                        float k,
                        int variant,
                        float xsize,
                        float ysize,
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
    fprintf(f, "mix= %d\n", mix);
    fclose(f);
}

/* Načte parametrický soubor .set. */
static void read_setup(const char *path,
                       int *n_cities,
                       unsigned long *max_cycles,
                       int *print_period,
                       float *t0,
                       float *k,
                       int *variant,
                       float *xsize,
                       float *ysize,
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
    if (fscanf(f, "mix= %d\n", mix) != 1) {
        fprintf(stderr, "Nelze přečíst mix\n");
        exit(EXIT_FAILURE);
    }

    fclose(f);
}

/* Načte konfiguraci ze souboru .conf. */
static void read_config(const char *path, int n_last, float *x, float *y) {
    if (!file_exists(path)) {
        fprintf(stderr, "Soubor neexistuje: %s\n", path);
        exit(EXIT_FAILURE);
    }

    FILE *f = xfopen(path, "r");
    int i = 0;
    while (i <= n_last && fscanf(f, "%f\t%f", &x[i], &y[i]) == 2) {
        i++;
    }
    fclose(f);

    if (i != n_last + 1) {
        fprintf(stderr, "Varování: načteno %d bodů, očekáváno %d\n", i, n_last + 1);
    }
}

/* Vytvoří textový .mol soubor (topologie řetězce pro vizualizaci). */
static void write_mol_header(const char *path,
                            int n_last,
                            unsigned long max_cycles,
                            int print_period,
                            float t0,
                            float k) {
    FILE *molf = xfopen(path, "a");
    const int n_atoms = n_last + 1;

    fprintf(molf, "! Soubor vygenerovan v ramci reseni problemu obchodniho cestujiciho.\n");
    fprintf(molf,
            "!Number of cities=%d, Steps=%lu, Printed every %d steps, Initial Temperature= %f, k=%f\n",
            n_atoms,
            max_cycles,
            print_period,
            t0,
            k);
    fprintf(molf, "\nnumber_of_atoms = %d\n\n\n", n_atoms);
    fprintf(molf, "atoms\n");
    fprintf(molf, "! i  atom-id  a-type charge chir nbonds bound_atoms\n");

    /* 0 je navázaná na 1 a N */
    fprintf(molf, "%3d %3d-H \t H      0    0      2     %d %d\n", 0, 0, 1, n_last);

    for (int i = 1; i < n_last; i++) {
        fprintf(molf, "%3d %3d-H \t H      0    0      2     %d %d\n", i, i, i - 1, i + 1);
    }

    /* N je navázaná na N-1 a 0 */
    fprintf(molf, "%3d %3d-H \t H      0    0      2     %d %d\n", n_last, n_last, n_last - 1, 0);

    fclose(molf);
}

/* Zapíše binární .plb hlavičku (2 floaty). */
static void write_plb_header(const char *path, const float header[2]) {
    FILE *plb = xfopen(path, "wb");
    fwrite(header, sizeof(float), 2, plb);
    fclose(plb);
}

/* Zapíše do .trip jeden záznam (krok, frame, teplota, délka). */
static void append_trip(const char *path,
                        float length,
                        unsigned long step,
                        int print_period,
                        float temperature) {
    FILE *f = xfopen(path, "a");
    fprintf(f, "%15lu\t%lu\t%f\t%f\n", step, (unsigned long)(step / (unsigned long)print_period),
            temperature, fabsf(length));
    fclose(f);
}

/* Vypíše dobu běhu do .trip a zároveň na stdout. */
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
/* Náhodná čísla                                                             */
/* ------------------------------------------------------------------------- */

/* Náhodné reálné číslo v [0,1]. */
static float rand01(void) {
    return (float)rand() / (float)RAND_MAX;
}

/* Náhodné reálné číslo v [a,b]. */
static float rand_range(float a, float b) {
    return rand01() * (b - a) + a;
}

/* Náhodné celé číslo v [0,n_last]. */
static int rand_index(int n_last) {
    return rand() % (n_last + 1);
}

/* Náhodně promíchá pořadí měst (Fisher–Yates). */
static void shuffle_points(int n_last, float *x, float *y) {
    /* n_last je index posledního prvku, počet prvků = n_last+1 */
    const int n = n_last + 1;

    int *perm = (int *)malloc((size_t)n * sizeof(int));
    float *x0 = (float *)malloc((size_t)n * sizeof(float));
    float *y0 = (float *)malloc((size_t)n * sizeof(float));
    if (!perm || !x0 || !y0) {
        fprintf(stderr, "Nedostatek paměti\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < n; i++) {
        perm[i] = i;
        x0[i] = x[i];
        y0[i] = y[i];
    }

    for (int i = n - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        int t = perm[i];
        perm[i] = perm[j];
        perm[j] = t;
    }

    for (int i = 0; i < n; i++) {
        x[i] = x0[perm[i]];
        y[i] = y0[perm[i]];
    }

    free(perm);
    free(x0);
    free(y0);
}

/* ------------------------------------------------------------------------- */
/* Geometrie a délka cesty                                                   */
/* ------------------------------------------------------------------------- */

/* Euklidovská vzdálenost dvou měst. */
static float dist_xy(const float *x, const float *y, int i, int j) {
    const float dx = x[i] - x[j];
    const float dy = y[i] - y[j];
    return sqrtf(dx * dx + dy * dy);
}

/* Délka hrany i -> i+1 (resp. poslední -> první). */
static float edge_length(const float *x, const float *y, int i, int n_last) {
    if (i == n_last) {
        return dist_xy(x, y, n_last, 0);
    }
    return dist_xy(x, y, i, i + 1);
}

/* Celková délka okružní cesty. */
static float tour_length(const float *x, const float *y, int n_last) {
    float sum = 0.0f;
    for (int i = 0; i <= n_last; i++) {
        sum += edge_length(x, y, i, n_last);
    }
    return sum;
}

/*
 * Změna délky cesty při otočení úseku [cb1, cb2].
 * Jde o standardní 2-opt:
 *   ... (cb1-1)--cb1 ... cb2--(cb2+1) ...
 * po otočení:
 *   ... (cb1-1)--cb2 ... cb1--(cb2+1) ...
 */
static float delta_2opt(const float *x, const float *y, int cb1, int cb2, int n_last) {
    const int a = cb1 - 1;
    const int b = cb1;
    const int c = cb2;
    const int d = cb2 + 1;

    const int aa = (a < 0) ? n_last : a;
    const int dd = (d > n_last) ? 0 : d;

    const float old_len = dist_xy(x, y, aa, b) + dist_xy(x, y, c, dd);
    const float new_len = dist_xy(x, y, aa, c) + dist_xy(x, y, b, dd);
    return new_len - old_len;
}

/* Provede in-place otočení úseku [cb1, cb2] v polích x,y. */
static void apply_2opt(float *x, float *y, int cb1, int cb2) {
    while (cb1 < cb2) {
        float tmp = x[cb1];
        x[cb1] = x[cb2];
        x[cb2] = tmp;

        tmp = y[cb1];
        y[cb1] = y[cb2];
        y[cb2] = tmp;

        cb1++;
        cb2--;
    }
}

/* ------------------------------------------------------------------------- */
/* "MojeOprava": lokální prohledání permutací v okně POCET                    */
/* ------------------------------------------------------------------------- */

static int factorial(int n) {
    if (n <= 1) {
        return 1;
    }
    return n * factorial(n - 1);
}

/*
 * Naplní pole a[] postupností POCET měst začínající od 'city' (s obtočením).
 * n_last je index posledního města.
 */
static void window_cities(int n_last, int city, int *a) {
    for (int i = 0; i < POCET; i++) {
        int v = city + i;
        if (v > n_last) {
            v -= (n_last + 1);
        }
        if (v < 0) {
            v += (n_last + 1);
        }
        a[i] = v;
    }
}

/*
 * Vrátí "kruhový" index: pokud i vyjede mimo [0,n_last], obtočí se.
 */
static int wrap_index(int i, int n_last) {
    if (i >= 0) {
        return (i <= n_last) ? i : (i - n_last - 1);
    }
    return n_last + i + 1;
}

/* Najde index hodnoty v poli délky POCET, pokud nenajde vrátí -1. */
static int find_in_window(const int *arr, int value) {
    for (int i = 0; i < POCET; i++) {
        if (arr[i] == value) {
            return i;
        }
    }
    return -1;
}

/*
 * Spočítá délku úseku: (city-1) -> okno POCET -> (city+POCET).
 * 'order' je permutace POCET měst v okně.
 */
static float window_path_length(const int *order,
                                const float *x,
                                const float *y,
                                int city,
                                int n_last) {
    int list[POCET + 2];
    list[0] = wrap_index(city - 1, n_last);
    list[POCET + 1] = wrap_index(city + POCET, n_last);

    for (int i = 0; i < POCET; i++) {
        list[i + 1] = order[i];
    }

    float sum = 0.0f;
    for (int i = 0; i < POCET + 1; i++) {
        sum += dist_xy(x, y, list[i], list[i + 1]);
    }
    return sum;
}

/*
 * Vygeneruje N-tou permutaci (lexikograficky) indexů 0..POCET-1 a aplikuje ji
 * na okno měst začínající v 'city'. Výsledek uloží do D[] (konkrétní čísla měst).
 */
static void nth_permutation_window(int n_last, int nth, int *D, int city) {
    const int fact = factorial(POCET);
    if (nth < 0 || nth >= fact) {
        fprintf(stderr,
                "nth_permutation_window: poradi mimo rozsah (0..%d), nth=%d\n",
                fact - 1,
                nth);
        nth %= fact;
        if (nth < 0) {
            nth += fact;
        }
    }

    int base[POCET];
    window_cities(n_last, city, base);

    /* a[] = permutace indexů 0..POCET-1 */
    int a[POCET];
    for (int i = 0; i < POCET; i++) {
        a[i] = i;
    }

    if (nth == 0) {
        for (int i = 0; i < POCET; i++) {
            D[i] = base[i];
        }
        return;
    }

    /* Lexikografické generování permutací až do pořadí nth. */
    for (int step = 0; step < nth; step++) {
        int i = POCET - 2;
        while (i >= 0 && a[i] > a[i + 1]) {
            i--;
        }
        int j = POCET - 1;
        while (a[j] < a[i]) {
            j--;
        }
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

    for (int i = 0; i < POCET; i++) {
        D[i] = base[a[i]];
    }
}

/*
 * Lokálně opraví pořadí měst v okně POCET začínajícím v 'city':
 *  - vyzkouší všechny permutace okna
 *  - vybere nejlepší podle varianty (min/max)
 *  - přepíše odpovídající segment v X,Y
 */
static void improve_city_window(float *x, float *y, int city, int n_last, int variant) {
    const int fact = factorial(POCET);

    int best_order[POCET];
    float best_value = 0.0f;

    for (int i = 0; i < fact; i++) {
        int order[POCET];
        nth_permutation_window(n_last, i, order, city);

        float val = (float)variant * window_path_length(order, x, y, city, n_last);
        if (i == 0 || val <= best_value) {
            best_value = val;
            memcpy(best_order, order, sizeof(best_order));
        }
    }

    /* Mapa z "aktuálních měst v okně" -> nové pořadí */
    int window_ids[POCET];
    for (int i = 0; i < POCET; i++) {
        window_ids[i] = wrap_index(city + i, n_last);
    }

    float x_tmp[POCET];
    float y_tmp[POCET];
    for (int i = 0; i < POCET; i++) {
        x_tmp[i] = x[window_ids[i]];
        y_tmp[i] = y[window_ids[i]];
    }

    for (int i = 0; i < POCET; i++) {
        const int idx_in_old = find_in_window(window_ids, best_order[i]);
        if (idx_in_old < 0) {
            continue; /* nemělo by nastat */
        }
        x[window_ids[i]] = x_tmp[idx_in_old];
        y[window_ids[i]] = y_tmp[idx_in_old];
    }
}

/* Provede lokální opravu pro všechna města a vrátí dobu běhu (s). */
static double my_correction(float *x, float *y, int n_last, int variant) {
    clock_t start = clock();
    for (int i = 0; i <= n_last; i++) {
        improve_city_window(x, y, i, n_last, variant);
    }
    clock_t finish = clock();
    return (double)(finish - start) / (double)CLOCKS_PER_SEC;
}


/* ------------------------------------------------------------------------- */
/* Simulované žíhání (SA)                                                    */
/* ------------------------------------------------------------------------- */

/*
 * Provede simulované žíhání nad aktuálním pořadím měst (X[], Y[]).
 *
 * Vstupy:
 *  - X, Y: pole souřadnic měst v pořadí trasy (cyklus 0..n_last..0)
 *  - n_last: poslední platný index (N-1)
 *  - variant: 1 pro minimalizaci délky, -1 pro maximalizaci (přenásobí energii)
 *  - max_cycles: počet cyklů
 *  - print_period: perioda zápisu snímků a logu
 *  - t0, k: parametry chladicího plánu
 *  - xsize, ysize: rozměry boxu (použito pro korekci pravděpodobnosti)
 *  - z: Z-ová souřadnice pro výstup do .plb
 *  - box_size: 3 floaty zapisované před každým framem do .plb
 *  - path_plb, path_trip: výstupní soubory (append)
 *
 * In/Out:
 *  - *temperature: na vstupu počáteční teplota, na výstupu poslední použitá
 *  - *cycle: na vstupu start (typicky 1), na výstupu hodnota po dokončení (max_cycles+1)
 *
 * Návrat:
 *  - doba běhu SA v sekundách
 */

static inline uint32_t xorshift32(uint32_t *state) {
    /* Velmi rychlý RNG; state nesmí být 0 */
    uint32_t x = *state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    *state = x;
    return x;
}

static inline float rng01(uint32_t *state) {
    /* Uniformní [0,1); 24 bitů mantisy je pro SA dostačující */
    return (xorshift32(state) >> 8) * (1.0f / 16777216.0f); /* 2^24 */
}

static inline int rng_index(uint32_t *state, int n) {
    /* Index 0..n-1 */
    return (int)(xorshift32(state) % (uint32_t)n);
}

static inline int imod(int a, int n) {
    a %= n;
    return (a < 0) ? a + n : a;
}

static inline float dist_at(const float *dist, int N, int a, int b) {
    return dist[a * N + b];
}

/* Δ délky při swapu prvků tour[i] a tour[j] v cyklické trase */
static inline float delta_swap(const int *tour, const float *dist, int N, int i, int j) {
    if (i == j) return 0.0f;
    if (j < i) {
        int t = i;
        i = j;
        j = t;
    }

    const int im1 = imod(i - 1, N);
    const int ip1 = imod(i + 1, N);
    const int jm1 = imod(j - 1, N);
    const int jp1 = imod(j + 1, N);

    const int A = tour[im1], B = tour[i], C = tour[ip1];
    const int D = tour[jm1], E = tour[j], F = tour[jp1];

    float oldc = 0.0f, newc = 0.0f;

    /* Speciální případ: i a j jsou sousedé (i -> j) */
    if (ip1 == j) {
        /* ... A - B - E - F ...  =>  ... A - E - B - F ... */
        oldc = dist_at(dist, N, A, B) + dist_at(dist, N, B, E) + dist_at(dist, N, E, F);
        newc = dist_at(dist, N, A, E) + dist_at(dist, N, E, B) + dist_at(dist, N, B, F);
        return newc - oldc;
    }

    /* Speciální případ: sousedé přes konec cyklu (j -> i) */
    if (jp1 == i) {
        /* ... D - E - B - C ...  =>  ... D - B - E - C ... */
        oldc = dist_at(dist, N, D, E) + dist_at(dist, N, E, B) + dist_at(dist, N, B, C);
        newc = dist_at(dist, N, D, B) + dist_at(dist, N, B, E) + dist_at(dist, N, E, C);
        return newc - oldc;
    }

    /* Obecný případ: mění se 4 hrany */
    oldc = dist_at(dist, N, A, B) + dist_at(dist, N, B, C) + dist_at(dist, N, D, E) + dist_at(dist, N, E, F);
    newc = dist_at(dist, N, A, E) + dist_at(dist, N, E, C) + dist_at(dist, N, D, B) + dist_at(dist, N, B, F);

    return newc - oldc;
}

static inline void apply_swap(int *tour, int i, int j) {
    const int t = tour[i];
    tour[i] = tour[j];
    tour[j] = t;
}

/*
 * Provede simulované žíhání nad aktuálním pořadím měst.
 *
 * ULTRAFAST (varianta A): návrh tahu = swap dvou měst (O(1) aplikace).
 * Klíčové optimalizace:
 *  - trasa je reprezentována permutací `tour[]`, souřadnice měst jsou fixní (`px[]`,`py[]`)
 *  - předpočtená matice vzdáleností `dist[N*N]` (float, 1D pole kvůli cache)
 *  - inkrementální délka trasy `current_len` + delta výpočet O(1)
 *  - expf nahrazen lookup tabulkou exp(-x) pro x∈[0, XMAX] (rychlejší než expf)
 *  - výstupní soubory jsou otevřené po dobu běhu SA (minimalizace I/O overhead)
 */
static double run_simulated_annealing(float *X,
                                      float *Y,
                                      int n_last,
                                      int variant,
                                      unsigned long max_cycles,
                                      int print_period,
                                      float t0,
                                      float k,
                                      float xsize,
                                      float ysize,
                                      float z,
                                      const float box_size[3],
                                      const char *path_plb,
                                      const char *path_trip,
                                      float *temperature,
                                      unsigned long *cycle) {
    const int N = n_last + 1;
    const float korr = sqrtf(xsize * ysize);

    /* Fixní souřadnice měst a počáteční permutace */
    float *px = (float *)malloc((size_t)N * sizeof(float));
    float *py = (float *)malloc((size_t)N * sizeof(float));
    int *tour = (int *)malloc((size_t)N * sizeof(int));

    if (!px || !py || !tour) {
        fprintf(stderr, "Memory allocation failed (N=%d)\n", N);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < N; i++) {
        px[i] = X[i];
        py[i] = Y[i];
        tour[i] = i;
    }

    /* Matice vzdáleností */
    float *dist = (float *)malloc((size_t)N * (size_t)N * sizeof(float));
    if (!dist) {
        fprintf(stderr, "Memory allocation failed for dist matrix (N=%d)\n", N);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < N; i++) {
        dist[(size_t)i * N + i] = 0.0f;
        for (int j = i + 1; j < N; j++) {
            const float dx = px[i] - px[j];
            const float dy = py[i] - py[j];
            const float d = sqrtf(dx * dx + dy * dy);
            dist[(size_t)i * N + j] = d;
            dist[(size_t)j * N + i] = d;
        }
    }

    /* Výpočet počáteční délky (cyklus) */
    float current_len = 0.0f;
    for (int i = 0; i < N - 1; i++) {
        current_len += dist_at(dist, N, tour[i], tour[i + 1]);
    }
    current_len += dist_at(dist, N, tour[N - 1], tour[0]);
    current_len *= (float)variant;

    /* Lookup tabulka pro exp(-x), x∈[0, XMAX] */
    enum { EXP_TABLE_SIZE = 65536 };
    const float XMAX = 20.0f;
    const float inv_step = (float)EXP_TABLE_SIZE / XMAX;

    float *expLUT = (float *)malloc(((size_t)EXP_TABLE_SIZE + 1) * sizeof(float));
    if (!expLUT) {
        fprintf(stderr, "Memory allocation failed for exp lookup table\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i <= EXP_TABLE_SIZE; i++) {
        const float x = (float)i / inv_step;
        expLUT[i] = expf(-x);
    }

    /* Otevři výstupy jednou */
    FILE *plb = xfopen(path_plb, "ab");
    FILE *trip = xfopen(path_trip, "a");

    /* RNG seed (nesmí být 0) */
    uint32_t rng_state = (uint32_t)time(NULL) ^ 0x9E3779B9u;
    if (rng_state == 0) rng_state = 1u;

    clock_t start = clock();

    while (*cycle <= max_cycles) {
        /* Chladicí plán */
        *temperature = max_cycles * k * t0 / (*cycle + (max_cycles * k - 1.0f));
        const float beta = 1.0f / ((*temperature) * korr); /* 1/(T*korr) */

        /* Návrh: swap dvou pozic v trase */
        int i = rng_index(&rng_state, N);
        int j = rng_index(&rng_state, N);
        if (i == j) {
            (*cycle)++;
            continue;
        }

        const float delta_len = (float)variant * delta_swap(tour, dist, N, i, j);

        /* Přijetí tahu */
        int accept = 0;
        if (delta_len <= 0.0f) {
            accept = 1;
        } else {
            const float x = delta_len * beta; /* dE/(T*korr) */
            if (x < XMAX) {
                /* LUT + lineární interpolace */
                const float fx = x * inv_step;
                const int idx = (int)fx;
                const float t = fx - (float)idx;
                const float p = expLUT[idx] + t * (expLUT[idx + 1] - expLUT[idx]);
                if (rng01(&rng_state) < p) accept = 1;
            } else {
                accept = 0; /* prakticky nulová pravděpodobnost */
            }
        }

        if (accept) {
            apply_swap(tour, i, j);
            current_len += delta_len;
        }

        if ((*cycle) % (unsigned long)print_period == 0) {
            /* .plb frame */
            fwrite(box_size, sizeof(float), 3, plb);
            for (int p = 0; p < N; p++) {
                const int city = tour[p];
                const float xx = px[city];
                const float yy = py[city];
                fwrite(&xx, sizeof(float), 1, plb);
                fwrite(&yy, sizeof(float), 1, plb);
                fwrite(&z, sizeof(float), 1, plb);
            }

            /* .trip log */
            fprintf(trip, "%f %lu %d %f\n", current_len, *cycle, print_period, *temperature);
        }

        (*cycle)++;
    }

    /* Zapiš finální frame pokud nebyl přesně na periodu */
    if (max_cycles % (unsigned long)print_period != 0) {
        fwrite(box_size, sizeof(float), 3, plb);
        for (int p = 0; p < N; p++) {
            const int city = tour[p];
            const float xx = px[city];
            const float yy = py[city];
            fwrite(&xx, sizeof(float), 1, plb);
            fwrite(&yy, sizeof(float), 1, plb);
            fwrite(&z, sizeof(float), 1, plb);
        }
        fprintf(trip, "%f %lu %d %f\n", current_len, max_cycles, print_period, *temperature);
    }

    fclose(plb);
    fclose(trip);

    /* Přepiš výsledek zpět do X,Y v pořadí trasy (kvůli dalším částem programu) */
    for (int p = 0; p < N; p++) {
        const int city = tour[p];
        X[p] = px[city];
        Y[p] = py[city];
    }

    free(expLUT);
    free(dist);
    free(tour);
    free(px);
    free(py);

    clock_t finish = clock();
    return (double)(finish - start) / (double)CLOCKS_PER_SEC;
}


/* ------------------------------------------------------------------------- */
/* Hlavní program                                                            */
/* ------------------------------------------------------------------------- */

static void print_usage(const char *prog) {
    printf("WRONG number of arguments.\n");
    printf("Program for solving Travelling salesman problem by simulated annealing.\n");
    printf("Number of cities must be higher than 4.\n\n");
    printf("Usage:\n");
    printf("  %s N NAME Cyc Pr T0 k xsize ysize [maximize]\n", prog);
    printf("  If optional last argument is present, distance will be MAXIMIZED.\n\n");
    printf("  N      number of cities\n");
    printf("  NAME   base name of output files\n");
    printf("  Cyc    number of cycles\n");
    printf("  Pr     print period in cycles\n");
    printf("  T0     initial temperature\n");
    printf("  k      small constant (e.g. 0.0001), larger => slower cooling\n");
    printf("  xsize, ysize  size of simulation box\n\n");
    printf("Example:\n");
    printf("  %s 50 Vystup 5000000 20000 200 0.0001 1 1\n\n", prog);
    printf("Alternative:\n");
    printf("  %s NAME\n", prog);
    printf("  Files NAME.conf and NAME.set must exist.\n\n");
    printf("WARNING:\n");
    printf("  NAME.mol NAME.plb and NAME.trip will be overwritten if exists!\n\n");
}

int main(int argc, char *argv[]) {
    /* N = počet měst, v kódu často pracujeme s n_last = N-1 (poslední index) */
    int N = 0;
    int print_period = 100;
    int variant = 1; /* 1 = minimalizace, -1 = maximalizace */
    int mix = 1;

    unsigned long max_cycles = 5000;
    unsigned long cycle = 1;

    float xsize = 1.0f, ysize = 1.0f;
    float z = 0.0f;
    float t0 = 0.0f, k = 0.0f;

    /* Náhodné semínko */
    srand((unsigned)time(NULL));

    if ((argc != 9) && (argc != 10) && (argc != 2)) {
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
    snprintf(path_mol, sizeof(path_mol), "trajectory.mol"); /* původně fixně trajectory.mol */
    snprintf(path_fin, sizeof(path_fin), "%s.fin", base);

    float temperature = 0.0f;

    /* Načtení/parsování argumentů */
    if (argc == 2) {
        read_setup(path_set, &N, &max_cycles, &print_period, &t0, &k, &variant, &xsize, &ysize, &mix);
    } else {
        sscanf(argv[1], "%d", &N);
        sscanf(argv[3], "%lu", &max_cycles);
        sscanf(argv[4], "%d", &print_period);
        sscanf(argv[5], "%f", &t0);
        sscanf(argv[6], "%f", &k);
        sscanf(argv[7], "%f", &xsize);
        sscanf(argv[8], "%f", &ysize);
        if (argc == 10) {
            variant = -1;
        }
    }

    if (N <= 4) {
        fprintf(stderr, "Number of cities must be higher than 4 (N=%d)\n", N);
        return EXIT_FAILURE;
    }

    /* n_last = poslední index v polích (0..n_last), počet prvků = N */
    const int n_last = N - 1;

    temperature = t0;

    /* Hlavička pro .plb */
    float box_size[3] = {xsize, ysize, z};
    float header[2] = {(float)N, -3.0f};

    /* Alokace (N prvků => index 0..N-1) */
    float *X = (float *)malloc((size_t)N * sizeof(float));
    float *Y = (float *)malloc((size_t)N * sizeof(float));
    if (!X || !Y) {
        fprintf(stderr, "Nedostatek paměti\n");
        return EXIT_FAILURE;
    }

    if (argc == 2) {
        /* Načti konfiguraci ze souboru */
        read_config(path_conf, n_last, X, Y);

        if (mix != 0) {
            shuffle_points(n_last, X, Y);
        }
    } else {
        /* Vygeneruj náhodnou konfiguraci */
        for (int i = 0; i <= n_last; i++) {
            X[i] = rand_range(0.0f, xsize);
            Y[i] = rand_range(0.0f, ysize);
        }
    }

    /* Připrav výstupy */
    write_config(path_conf, n_last, X, Y);
    write_setup(path_set, N, max_cycles, print_period, t0, k, variant, xsize, ysize, mix);

    write_plb_header(path_plb, header);
    write_mol_header(path_mol, n_last, max_cycles, print_period, t0, k);

    /* Inicializuj .trip */
    {
        FILE *trip = xfopen(path_trip, "w");
        fprintf(trip,
                "!Number of cities=%d, Steps=%lu, Printed every %d steps, Initial Temperature= %f, k=%f\n",
                N,
                max_cycles,
                print_period,
                t0,
                k);
        fprintf(trip, "!Steps\tFrames\tTemperature\tLength of trip\n");
        fclose(trip);
    }

    /* Zapiš počáteční snímek do .plb */
    {
        FILE *plb = xfopen(path_plb, "ab");
        fwrite(box_size, sizeof(float), 3, plb);
        for (int i = 0; i <= n_last; i++) {
            fwrite(&X[i], sizeof(float), 1, plb);
            fwrite(&Y[i], sizeof(float), 1, plb);
            fwrite(&z, sizeof(float), 1, plb);
        }
        fclose(plb);
    }

    append_trip(path_trip, (float)variant * tour_length(X, Y, n_last), cycle, print_period, temperature);

    /* Korekce pro škálování pravděpodobnosti (původně sqrt(x*y)) */

    /* Simulované žíhání */
    const double sa_time = run_simulated_annealing(X,
                                                   Y,
                                                   n_last,
                                                   variant,
                                                   max_cycles,
                                                   print_period,
                                                   t0,
                                                   k,
                                                   xsize,
                                                   ysize,
                                                   z,
                                                   box_size,
                                                   path_plb,
                                                   path_trip,
                                                   &temperature,
                                                   &cycle);

    append_runtime(path_trip, "Simulated annealing", sa_time);
    printf("Final trip length from simulated annealing is: %f\n", tour_length(X, Y, n_last));

    /* MojeOprava: několikrát zopakuj lokální opravu */
    double corr_time_total = 0.0;
    float final_sum = 0.0f;

    for (int r = 0; r <= 3; r++) {
        corr_time_total += my_correction(X, Y, n_last, variant);
        final_sum = tour_length(X, Y, n_last);

        FILE *plb = xfopen(path_plb, "ab");
        fwrite(box_size, sizeof(float), 3, plb);
        for (int i = 0; i <= n_last; i++) {
            fwrite(&X[i], sizeof(float), 1, plb);
            fwrite(&Y[i], sizeof(float), 1, plb);
            fwrite(&z, sizeof(float), 1, plb);
        }
        fclose(plb);

        append_trip(path_trip, final_sum, cycle, print_period, temperature);
    }

    append_runtime(path_trip, "My correction", corr_time_total);
    {
        FILE *trip = xfopen(path_trip, "a");
        fprintf(trip, "!Trip length after my correction follows:\n");
        fclose(trip);
    }
    printf("Final length after My correction is: %f\n", final_sum);

    write_config(path_fin, n_last, X, Y);

    free(X);
    free(Y);
    return EXIT_SUCCESS;
}
