[1mdiff --git a/src/mmgroup/dev/generators/gen_ufind_lin2.ske b/src/mmgroup/dev/generators/gen_ufind_lin2.ske[m
[1mindex 214e0d0..fd50088 100644[m
[1m--- a/src/mmgroup/dev/generators/gen_ufind_lin2.ske[m
[1m+++ b/src/mmgroup/dev/generators/gen_ufind_lin2.ske[m
[36m@@ -77,7 +77,7 @@[m [mcheck_out_buf32(uint32_t *p_buf, uint32_t len, uint32_t min_len)[m
 // A structure for the data in the array ``a`` in function[m
 // ``gen_ufind_lin2_init``[m
 typedef struct {[m
[31m-   uint32_t status;    // Status of the array[m
[32m+[m[32m   int32_t status;     // Status of the array[m
    uint32_t n;         // dimenstion ``n`` of ``GF(2)^n``[m
    uint32_t n_max_g;   // Maximum No of generators of the group ``G``[m
    uint32_t n_g;       // Number of generators of the group ``G``[m
[36m@@ -93,9 +93,10 @@[m [mtypedef struct {[m
 [m
 [m
 [m
[31m-static inline void load_lin2_info(uint32_t *a, lin2_type *ps)[m
[32m+[m[32mstatic inline int32_t load_lin2_info(uint32_t *a, lin2_type *ps)[m
 {[m
[31m-    ps->status = a[0];[m
[32m+[m[32m    if (a == NULL) return ERR_GEN_UFIND_INT_LIN2 - 1;[m
[32m+[m[32m    ps->status = (int32_t)a[0];[m
     ps->n = a[1];[m
     ps->n_max_g = a[2];[m
     ps->n_g = a[3];[m
[36m@@ -103,17 +104,19 @@[m [mstatic inline void load_lin2_info(uint32_t *a, lin2_type *ps)[m
     ps->p_t = a + 5;[m
     ps->p_o = ps->p_t + ((size_t)1UL << ps->n);[m
     ps->p_g = ps->p_o + ((size_t)1UL << ps->n);[m
[32m+[m[32m    return ps->status;[m
 }[m
 [m
 [m
 static inline int32_t store_lin2_info(lin2_type *ps, uint32_t *a)[m
 {[m
[32m+[m[32m    if ((int32_t)a[0] < 0) return (int32_t)a[0];[m
     a[0] = ps->status;[m
     a[1] = ps->n;[m
     a[2] = ps->n_max_g;[m
     a[3] = ps->n_g;[m
     a[4] = ps->n_orbits;[m
[31m-    return (int32_t) a[0];[m
[32m+[m[32m    return ps->status;[m
 }[m
 [m
 [m
[36m@@ -122,6 +125,7 @@[m [mstatic inline int32_t lin2_error(uint32_t *a, int32_t status)[m
 {[m
     int32_t status_a;[m
     if (status >= 0) return status;[m
[32m+[m[32m    if (a == NULL) return ERR_GEN_UFIND_INT_LIN2 - 1;[m
     status_a = (int32_t)(a[0]);[m
     if (status_a >= 0) a[0] = (uint32_t)(status);[m
     return status;[m
[36m@@ -241,9 +245,9 @@[m [mmat_inverse(uint32_t *m, uint32_t n, uint32_t *m_inv)[m
 // Store inverse of the ``n`` times ``n`` bit matrix ``m`` in the[m
 // array ``m_inv``. Return 0 if inverse can be computed and -1 if not.[m
 {[m
[31m-    uint64_t a[24];[m
[32m+[m[32m    uint64_t a[LIN2_MAX_N];[m
     uint32_t i, mask = (1UL << n) - 1;[m
[31m-    if (n == 0 || n > sizeof(a) / sizeof(uint64_t)) return -1;[m
[32m+[m[32m    if (n == 0 || n > LIN2_MAX_N) return ERR_GEN_UFIND_INVERSE;[m
     for (i = 0; i < n; ++i) a[i] = m[i];[m
     if (bitmatrix64_inv(a, n)) return ERR_GEN_UFIND_INVERSE;[m
     for (i = 0; i < n; ++i) m_inv[i] = (uint32_t)(a[i] & mask);[m
[36m@@ -322,21 +326,16 @@[m [mdone:[m
 [m
 [m
 [m
[31m-static inline int32_t finalize_initalization(uint32_t *a)[m
[32m+[m[32mstatic inline[m[41m [m
[32m+[m[32mint32_t finalize_initalization(uint32_t *a, lin2_type *ps)[m
 {[m
     lin2_type s;[m
[31m-    int32_t status;[m
[32m+[m[32m    int32_t status = load_lin2_info(a, ps);[m
     uint32_t i, mask, *p_ind = NULL;[m
 [m
[31m-    if (a == NULL) return  ERR_GEN_UFIND_INT_LIN2 - 11;[m
[31m-    if (a[0] > 0 && a[0] <= LIN2_MAX_STATUS) return a[0];[m
[31m-    [m
[32m+[m[32m    if (status != 0) return status;[m
     load_lin2_info(a, &s);[m
[31m-    if (a[0] > LIN2_MAX_STATUS) {[m
[31m-        if ((int32_t)(a[0]) < 0) return (int32_t)(a[0]);[m
[31m-        status = ERR_GEN_UFIND_INT_LIN2 - 12;[m
[31m-        goto done;[m
[31m-    }[m
[32m+[m
     // Finalize union-find algorithm on orbits in main table s.p_t[m
     status = gen_ufind_find_all_min(s.p_t, 1UL << s.n);[m
     if  (status < 1) {[m
[36m@@ -383,7 +382,7 @@[m [mstatic inline int32_t finalize_initalization(uint32_t *a)[m
     [m
 done:[m
     if (p_ind) free(p_ind);[m
[31m-    s.status = (uint32_t)(status);[m
[32m+[m[32m    ps->status = s.status = (uint32_t)(status);[m
     store_lin2_info(&s, a);[m
     return status;[m
 }[m
[36m@@ -439,8 +438,8 @@[m [min case of error.[m
 int32_t gen_ufind_lin2_init(uint32_t *a, uint32_t l_a, uint32_t n, uint32_t *g, uint32_t k)[m
 {[m
     int32_t l_a_expected, status;[m
[31m-    uint32_t i, *p_ind = NULL;[m
     lin2_type s;[m
[32m+[m[32m    uint32_t i;[m
 [m
     status = l_a_expected = gen_ufind_lin2_size(n, k);[m
     if (status < 0) goto done;[m
[36m@@ -451,21 +450,23 @@[m [mint32_t gen_ufind_lin2_init(uint32_t *a, uint32_t l_a, uint32_t n, uint32_t *g,[m
     s.n_max_g = k;[m
     s.n_g = 0;[m
     s.n_orbits = 0;[m
[31m-    store_lin2_info(&s, a);[m
[31m-    load_lin2_info(a, &s);[m
 [m
     // Perform union-find algorithm on orbits in main table s.p_t[m
     status = gen_ufind_init(s.p_t, 1UL << n);[m
     if (status < 0) goto done;[m
     store_lin2_info(&s, a);[m
[32m+[m
[32m+[m[32m    status = -333;[m
[32m+[m[32m    goto done;[m
[32m+[m
[32m+[m[32m    // Perform union-find algorithm on orbits in main table s.p_t[m
     for (i = 0; i < k; ++i) {[m
         status = load_bitmatrix(a, g + i * s.n, s.n, 1);[m
[31m-        if (status < 0) {status -= 1000; goto done;}[m
[32m+[m[32m        if (status < 0) goto done;[m
     }[m
[31m-    status = finalize_initalization(a);[m
[32m+[m[32m    status = finalize_initalization(a, &s);[m
 [m
 done:[m
[31m-    if (p_ind) free(p_ind);[m
     return lin2_error(a, status);[m
 }[m
 [m
[36m@@ -475,22 +476,6 @@[m [mdone:[m
 * Obtaining information from array dealing with a subgroup of SL_2(n)[m
 ************************************************************************/[m
 [m
[31m-/// @cond DO_NOT_DOCUMENT[m
[31m-[m
[31m-static inline int32_t[m
[31m-check_a_fast(uint32_t *a)[m
[31m-{[m
[31m-    if (a == NULL) return  ERR_GEN_UFIND_INT_LIN2 - 11;[m
[31m-    if (a[0] > LIN2_MAX_STATUS || a[0] == 0) {[m
[31m-        if ((int32_t)a[0] < 0) return (int32_t)a[0];[m
[31m-        a[0] = (uint32_t)(ERR_GEN_UFIND_INT_LIN2 - 12);[m
[31m-        return (int32_t)a[0];[m
[31m-    }[m
[31m-    return 0;[m
[31m-}[m
[31m-[m
[31m-/// @endcond[m
[31m-[m
 [m
 [m
 /** @brief Given a group acting on ``GF(2)^n`` the function returns ``n``[m
[36m@@ -502,8 +487,9 @@[m [mreturns the dimension of the vector space.[m
 // %%EXPORT px[m
 int32_t gen_ufind_lin2_dim(uint32_t *a)[m
 {[m
[31m-    int32_t status = check_a_fast(a);[m
[31m-    return status < 0 ? status : a[1];[m
[32m+[m[32m    lin2_type s;[m
[32m+[m[32m    int32_t status = load_lin2_info(a, &s);[m
[32m+[m[32m    return status < 0 ? status : s.n;[m
 }[m
 [m
 [m
[36m@@ -516,8 +502,9 @@[m [mreturns the number of generators of the group.[m
 // %%EXPORT px[m
 int32_t gen_ufind_lin2_n_gen(uint32_t *a)[m
 {[m
[31m-    int32_t status = check_a_fast(a);[m
[31m-    return status < 0 ? status : a[3];[m
[32m+[m[32m    lin2_type s;[m
[32m+[m[32m    int32_t status = load_lin2_info(a, &s);[m
[32m+[m[32m    return status < 0 ? status : s.n_g;[m
 }[m
 [m
 [m
[36m@@ -531,8 +518,9 @@[m [mof the group.[m
 // %%EXPORT px[m
 int32_t gen_ufind_lin2_n_orbits(uint32_t *a)[m
 {[m
[31m-    int32_t status = check_a_fast(a);[m
[31m-    return status < 0 ? status : a[4];[m
[32m+[m[32m    lin2_type s;[m
[32m+[m[32m    int32_t status = finalize_initalization(a, &s);[m
[32m+[m[32m    return status < 0 ? status : s.n_orbits;[m
 }[m
 [m
 [m
[36m@@ -556,9 +544,8 @@[m [min case of error; e.g. if the array ``g`` is too short.[m
 int32_t gen_ufind_lin2_gen(uint32_t *a, uint32_t i, uint32_t *g, uint32_t l_g)[m
 {[m
     lin2_type s;[m
[31m-    int32_t status = check_a_fast(a);[m
[32m+[m[32m    int32_t status = load_lin2_info(a, &s);[m
     if (status < 0) return status;[m
[31m-    load_lin2_info(a, &s);[m
     if (i >= 2 * s.n_g) return ERR_GEN_UFIND_IN_LARGE;[m
     if (check_out_buf32(g, l_g, s.n)) return ERR_GEN_UFIND_OUT_SHORT;[m
     memcpy(g, s.p_g + i * s.n, s.n * sizeof(uint32_t));[m
[36m@@ -581,14 +568,15 @@[m [mmay not be used![m
 int32_t gen_ufind_lin2_check(uint32_t *a, uint32_t len_a)[m
 {[m
     lin2_type s;[m
[31m-    int32_t a_size, status = check_a_fast(a);[m
[31m-    if (status < 0) return status;[m
[32m+[m[32m    int32_t a_size, status;[m
     if (len_a < 6) return lin2_error(a, ERR_GEN_UFIND_INT_LIN2 - 13);[m
[31m-    load_lin2_info(a, &s);[m
[32m+[m[32m    status = load_lin2_info(a, &s);[m
[32m+[m[32m    if (status < 0) return status;[m
     a_size = status = gen_ufind_lin2_size(s.n, s.n_g);[m
     if (status < 0) return lin2_error(a, status);[m
[31m-    if ((uint32_t)a_size > len_a) return lin2_error(a, ERR_GEN_UFIND_INT_LIN2 - 14);[m
[31m-    return a[0];[m
[32m+[m[32m    if ((uint32_t)a_size > len_a) return[m[41m [m
[32m+[m[32m        lin2_error(a, ERR_GEN_UFIND_INT_LIN2 - 14);[m
[32m+[m[32m    return status;[m
 }[m
 [m
 [m
[36m@@ -615,10 +603,9 @@[m [mint32_t gen_ufind_lin2_rep_v(uint32_t *a, uint32_t v)[m
 {[m
     lin2_type s;[m
     uint32_t mask, entry;[m
[31m-    int32_t status = check_a_fast(a);[m
[32m+[m[32m    int32_t status = finalize_initalization(a, &s);[m
 [m
     if (status < 0) return status;[m
[31m-    load_lin2_info(a, &s);[m
     mask = (1UL << s.n) - 1;[m
     v &= mask;[m
     entry = s.p_t[v];[m
[36m@@ -642,10 +629,9 @@[m [mint32_t gen_ufind_lin2_len_orbit_v(uint32_t *a, uint32_t v)[m
 {[m
     lin2_type s;[m
     uint32_t mask, rep, index_o;[m
[31m-    int32_t status = check_a_fast(a);[m
[32m+[m[32m    int32_t status = finalize_initalization(a, &s);[m
 [m
     if (status < 0) return status;[m
[31m-    load_lin2_info(a, &s);[m
     mask = (1UL << s.n) - 1;[m
     rep = gen_ufind_lin2_rep_v(a, v);[m
     index_o = s.p_t[rep] & mask;[m
[36m@@ -673,10 +659,9 @@[m [mint32_t gen_ufind_lin2_orbit_v(uint32_t *a, uint32_t v, uint32_t *r, uint32_t l_[m
     lin2_type s;[m
     uint32_t mask, rep, index_o, *p_o;[m
     int32_t length, i;[m
[31m-    int32_t status = check_a_fast(a);[m
[32m+[m[32m    int32_t status = finalize_initalization(a, &s);[m
 [m
     if (status < 0) return status;[m
[31m-    load_lin2_info(a, &s);[m
     mask = (1UL << s.n) - 1;[m
     rep = status = gen_ufind_lin2_rep_v(a, v);[m
     if (status < 0) return status;[m
[36m@@ -711,10 +696,9 @@[m [mint32_t gen_ufind_lin2_representatives(uint32_t *a, uint32_t *r, uint32_t l_r)[m
 {[m
     lin2_type s;[m
     uint32_t mask, index_o = 0, l_r1 = l_r, d;[m
[31m-    int32_t status = check_a_fast(a);[m
[32m+[m[32m    int32_t status = finalize_initalization(a, &s);[m
 [m
     if (status < 0) return status;[m
[31m-    load_lin2_info(a, &s);[m
     mask = (1UL << s.n) - 1;[m
 [m
     while (index_o < (1UL << s.n)) {[m
[36m@@ -759,16 +743,19 @@[m [mstore_gen(uint32_t *g, uint32_t n, uint32_t start, uint32_t maxlen, uint32_t *o)[m
 /** @brief Yet to be tested and documented![m
 [m
 */[m
[31m-static inline int32_t add_generators(uint32_t *a)[m
[32m+[m[32mstatic inline int32_t add_generators(uint32_t *a, lin2_type *p_s)[m
 {[m
     uint32_t *p_alloc = NULL;[m
     uint32_t *p_start, *p_end, *p_overflow, *p_g, *p_s_g, *p_g1, *p_b;[m
     uint32_t i, j, l_b;[m
[31m-    int32_t status = gen_ufind_lin2_check(a, 0x7ffffffUL);[m
     lin2_type s;[m
[32m+[m[32m    int32_t status = finalize_initalization(a, p_s);[m
 [m
     if (status < 0) return status;[m
[31m-    if (a[0] >= 2) return 0;[m
[32m+[m[32m    if (a[0] >= 2) {[m
[32m+[m[32m        a[0] = 2;[m
[32m+[m[32m        return 0;[m
[32m+[m[32m    }[m
     load_lin2_info(a, &s);[m
 [m
     // Allocate work buffer and set pointer to queue in work buffer[m
[36m@@ -829,12 +816,11 @@[m [mstatic inline int32_t add_generators(uint32_t *a)[m
             *p_end++ = img;[m
         }[m
     }[m
[31m-    a[0] = 2;  // update main status[m
[31m-    status = 0;[m
[32m+[m[32m    status = a[0] = p_s->status = 2;  // update main status[m
 done:[m
     if (p_alloc) free(p_alloc);[m
     if (status < 0) return lin2_error(a, status);[m
[31m-    return 0;[m
[32m+[m[32m    return status;[m
 }[m
 [m
 [m
[36m@@ -854,10 +840,9 @@[m [mDetails are yet to be documented![m
 int32_t gen_ufind_lin2_map_v_gen(uint32_t *a, uint32_t v)[m
 {[m
     uint32_t mask;[m
[31m-    int32_t status = add_generators(a);[m
     lin2_type s;[m
[32m+[m[32m    int32_t status = add_generators(a, &s);[m
     if (status < 0) return status;[m
[31m-    load_lin2_info(a, &s);[m
     mask = (1UL << s.n) - 1;[m
     return (s.p_t[v & mask] >> 24) & 0xff;[m
 }[m
[36m@@ -874,11 +859,10 @@[m [mDetails are yet to be documented![m
 // %%EXPORT px[m
 int32_t gen_ufind_lin2_map_v(uint32_t *a, uint32_t v, uint8_t *b, uint32_t l_b)[m
 {[m
[31m-    lin2_type s;[m
     uint32_t mask, w, l_b1, entry;[m
[31m-    int32_t status = add_generators(a);[m
[32m+[m[32m    lin2_type s;[m
[32m+[m[32m    int32_t status = add_generators(a, &s);[m
     if (status < 0) return status;[m
[31m-    load_lin2_info(a, &s);[m
     mask = (1UL << s.n) - 1;[m
     v &= mask;[m
     entry = s.p_t[v];[m
[36m@@ -936,10 +920,8 @@[m [mint32_t gen_ufind_lin2_finalize(uint32_t *a)[m
     lin2_type s;[m
     int32_t status, a_length, i;[m
     uint32_t a_sum = 0;[m
[31m-    if ((status = add_generators(a)) != 0) return status;[m
[31m-    load_lin2_info(a, &s);[m
[31m-    if (s.status != LIN2_MAX_STATUS) [m
[31m-        return status < 0 ? status :  ERR_GEN_UFIND_INT_LIN2 - 41;[m
[32m+[m[32m    status = add_generators(a, &s);[m
[32m+[m[32m    if (status < 0) return status;[m
     a_length =  gen_ufind_lin2_size(s.n, s.n_g);[m
     a[a_length - 2] = magic_content(&s);[m
     for (i = 0; i < a_length - 1; ++i) a_sum += a[i];[m
[36m@@ -962,10 +944,11 @@[m [mint32_t gen_ufind_lin2_check_finalized(uint32_t *a, uint32_t len_a)[m
     lin2_type s;[m
     int32_t status, a_length, i;[m
     uint32_t a_sum = 0;[m
[32m+[m[32m    status = load_lin2_info(a, &s);[m
[32m+[m[32m    if (status < 0) return status;[m
[32m+[m
     if ((status = gen_ufind_lin2_check(a, len_a)) < 0) return status;[m
[31m-    load_lin2_info(a, &s);[m
[31m-    if (s.status != LIN2_MAX_STATUS)[m
[31m-        return status < 0 ? status : ERR_GEN_UFIND_INT_LIN2 - 43;[m
[32m+[m[32m    if (status < LIN2_MAX_STATUS) return ERR_GEN_UFIND_INT_LIN2 - 43;[m
     a_length =  gen_ufind_lin2_size(s.n, s.n_g);[m
     if (a_length <= 0) return ERR_GEN_UFIND_INT_LIN2 - 44;;[m
     if (a[a_length - 2] != magic_content(&s)) return ERR_GEN_UFIND_INT_LIN2 - 45;[m
[36m@@ -1002,7 +985,7 @@[m [mint32_t gen_ufind_lin2_get_map(uint32_t *a, uint32_t *map, uint32_t l_map)[m
 {[m
     lin2_type s;[m
     uint32_t i;[m
[31m-    int32_t status = check_a_fast(a);[m
[32m+[m[32m    int32_t status = finalize_initalization(a, &s);[m
 [m
     if (status < 0) return status;[m
     load_lin2_info(a, &s);[m
[36m@@ -1030,10 +1013,9 @@[m [min case of failure, e.g. if the array ``t`` is too short.[m
 int32_t gen_ufind_lin2_get_table(uint32_t *a, uint32_t *t, uint32_t l_t)[m
 {[m
     lin2_type s;[m
[31m-    int32_t status = check_a_fast(a);[m
[32m+[m[32m    int32_t status = finalize_initalization(a, &s);[m
 [m
     if (status < 0) return status;[m
[31m-    load_lin2_info(a, &s);[m
     if (check_out_buf32(t, l_t, 1UL << s.n)) return ERR_GEN_UFIND_OUT_SHORT;[m
     memcpy(t, s.p_t, sizeof(uint32_t) << s.n);[m
     return 1UL << s.n;[m
[36m@@ -1067,10 +1049,9 @@[m [mint32_t gen_ufind_lin2_orbits(uint32_t *a, uint32_t *t, uint32_t l_t, uint32_t *[m
 {[m
     lin2_type s;[m
     uint32_t i;[m
[31m-    int32_t status = check_a_fast(a);[m
[32m+[m[32m    int32_t status = finalize_initalization(a, &s);[m
 [m
     if (status < 0) return status;[m
[31m-    load_lin2_info(a, &s);[m
     if (check_out_buf32(t, l_t, 1UL << s.n)) return ERR_GEN_UFIND_OUT_SHORT;[m
     if (check_out_buf32(x, l_x, s.n_orbits + 1)) return ERR_GEN_UFIND_OUT_SHORT;[m
     /* Ye olde waye to do the worke[m
[1mdiff --git a/src/mmgroup/tests/test_gen_xi/test_gen_ufind.py b/src/mmgroup/tests/test_gen_xi/test_gen_ufind.py[m
[1mindex b4f6c1e..ef93438 100644[m
[1m--- a/src/mmgroup/tests/test_gen_xi/test_gen_ufind.py[m
[1m+++ b/src/mmgroup/tests/test_gen_xi/test_gen_ufind.py[m
[36m@@ -194,6 +194,7 @@[m [mdef union_linear_high_level(generators):[m
     global a[m
     a = np.zeros(len_a, dtype = np.uint32)[m
     chk(gen_ufind_lin2_init(a, len_a, dim, gen.ravel(), n_gen))[m
[32m+[m[32m    1/0[m
     print("iii", a[:6], len(a), n_gen)[m
     t_len = 1 << chk(gen_ufind_lin2_dim(a))[m
     n_orbits = chk(gen_ufind_lin2_n_orbits(a))[m
