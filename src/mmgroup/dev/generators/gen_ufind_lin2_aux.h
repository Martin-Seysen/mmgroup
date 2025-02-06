/** @file gen_ufind_lin2_aux.h

Internal header file for files ``gen_ufind_lin2*.c``.

This file is used in the mmgroup poject 'as is', without preprocessing
by the code generator. The infomation stored in this file is not
relevant for the public C interface of the mmgroup project, so we
we do do export any documentation from this file.

The functions in file ``gen_ufind_lin2.c`` deal with the operation
of a group \f$G\f$ as a permutation group on the vector
space \f$V = \mbox{GF}_2^n\f$ for \f$n \leq 24\f$. The information
about the group \f$G\f$ is stored in an opaque array ``a``. The
internal structure of array ``a`` is documented in this headar file.
*/


#ifndef GEN_UFIND_LIN2_AUX_H
#define GEN_UFIND_LIN2_AUX_H


/************************************************************************
*  Stuctures for orbit arrays
************************************************************************/

/// @cond DO_NOT_DOCUMENT

/*********************************************************************
A function in module ``gen_ufind_lin2.c`` deals with a group \f$G\f$,
acting on a linear or affine space \f$V = GF_2^n\f$ and stores the
relevant information in an array ``a``, which is usually the first
parameter of such a function. This makes the interface to the outside
world simple; but it slightly complicates the internal operation of a
function in module ``gen_ufind_lin2.c``.

An array ``a`` must be generated with function ``gen_ufind_lin2_init``.
Any function dealing with that array should first call the (fast)
inline function ``load_lin2_info`` defined in this module. This function
copies header information and pointers to data from the array ``a`` into
a structure of type ``lin2_type``, for simplifying access to these data.
If that function changes any information in that structure then it
should copy the structure of type ``lin2_type`` back to the array ``a``
by calling function ``store_lin2_info``.

A structure ``s`` of type ``lin2_type`` contains an entry ``status``
which is  interpreted as follows.

  0: Initial status, generators of the group may be added.
  1: Deprecated, and no longer in use.
  2: A Schreier vector has been computed
     E.g. function ``gen_ufind_lin2_finalize`` upgrades to status 2
  LIN2_COMPRESSED: Array ``a`` has been compressed
     E.g. function ``gen_ufind_lin2_compress`` upgrades to this status
  Any negative status indicates an error.
     If an error is detected then the status will be set to an
     negative value. A function detecting an erroneous array ``a``
     must terminate with an error indication immediately.

Structure ``s`` has an entry ``n`` which is the dimension of \f$V\f$.

Structure ``s`` contains an entry ``p_t`` pointing to a table ``t``.

  In case ``status == 2`` this table has length ``1 << s.n``,
  where ``s.n`` is the dimension of the vector space \f$V\f$.
  Entry ``i`` in that table corresponds to vector ``i`` in \f$V\f$.
  The union-find algorithm in module ``gen_union_find.c`` is executed
  on this table in order to unite each vector in \f$V\f$ with its
  images under the operation of each generator of the group \f$G\f$.

  Bits 0,...,23 of entry ``i`` of array ``t`` is an index ``j``
  referring to table ``map`` referred by ``s.p_o``. Entry ``map[j]``
  describes the orbit of entry ``i``, as explained in the description
  of entry ``p_o``. This is valid if the status is 2.

  Bits 24,...,31 of entry ``i`` of array ``t`` is the number ``k`` of
  a generator ``g`` of the group \f$G\f$. Here a Schreier vector is a
  mapping of the vectors in  \f$V\f$ to the union of set of the
  generators of the group \f$G\f$ with the neutral element of \f$G\f$.
  The neutral element of \f$G\f$ is encoded as ``k = 0xfe``.

Structure ``s`` contains an entry ``p_o`` pointing to a table ``map``.
 
  In case ``status == 2`` this table has length ``1 << s.n``. It
  contains a partition of the vectors of \f$V\f$ into orbits under the
  action of the group \f$G\f$. The entries of each orbit are ordered.
  All orbits are stored contiguously in the table ``map``, ordered by
  the least vector in the orbit. Bits 0,...,23 of an entry describe the
  vector in \f$V\f$ stored in that entry. Bits 24,...,31 of the first
  entry of any orbit encode the length of that orbit, as given by
  function ``write_length_info``.
 
Structure ``s`` has an entry ``p_g`` pointing to generators of \f$G\f$.

  This array has ``lin2_generator_size(p_g->n, p_g->n_max_g)`` entries.
  Here ``p_g->n`` is the dimension of the space \f$V\f$,
  and ``p_g->n_max_g`` is the maximum number of generators that may be
  stored in the array.  ``p_g->n_g`` is the number of generators
  actually stored in the array.

  Function ``lin2_generator(&s, 2*i)`` returns a pointer to
  the ``i`` -th   generator. Function ``lin2_generator(&s, 2*i + 1)``
  returns a pointer to the inverse of that generator.

  Generators are stored in the order as they are entered (and accepted)
  by function ``gen_ufind_lin2_add``.

Entries ``s.p_t`` and ``s.p_o`` in case ``s.status = LIN2_COMPRESSED``

  If ``s.status`` is LIN2_COMPRESSED then entries ``s.p_t``
  and ``s.p_o`` contain different data as described in the sequel.

  In this case a subset of \f$V\f$ equal to a union of orbits is
  stored in the corresponding tables.

  Entry ``s.p_t`` points to a table ``t`` of length ``ps->n_vectors``
  containing a sorted list of all elements of \f$V\f$ stored in the
  compressed array. An entry in that table stores the corresponding
  element of \f$V\f$ in bits 8,...,31. Bits 0,...,7 of that entry is
  the number ``k`` of a generator ``g`` of the group \f$G\f$; these
  bits encode a Schreier vector as descibed above.

  Entry ``s.p_o`` contains a table of length ``2 * ps->n_orbits``,
  where ``ps->n_orbits`` is the the number of orbits of \f$G\f$
  on \f$V\f$ contained in the compressed array. Entry ``i`` of
  that table is an index of an entry in the table ``t`` containing
  the first element of that orbit.  Entry ``i + ps->n_orbits`` of
  that table is the length of that orbit.

*********************************************************************/

// A structure for interpreting the data in the array ``a`` in a
// function in module .
typedef struct {
   int32_t status;     // Status of the array
   uint32_t n;         // dimension ``n`` of ``GF(2)^n``
   uint32_t n_max_g;   // Maximum No of generators of the group ``G``
   uint32_t n_g;       // Number of generators of the group ``G``
   uint32_t n_orbits;  // No of ``G`` orbits in ``GF(2)^n``
   uint32_t n_vectors; // No of vectors, usually 1 << n
   uint32_t *p_t;      // Pointer to main table ``t`` inside ``a``
   uint32_t *p_o;      // Pointer to main table ``map`` inside ``a``
   uint32_t *p_g;      // Pointer to list of generators of group in ``a``
} lin2_type;

#define LIN2_LEN_HEADER  6UL  // Length of header of structrue above
#define LIN2_MAX_STATUS  2UL  // Maximum good status
#define LIN2_MAX_N      24UL  // Maximum dimension n of GF(2)^n
#define LIN2_MAX_N_G   127UL  // Maximum No n_g of generators of group

#define LIN2_COMPRESSED 0x10UL  // Status indicating a compressed orbit array



/************************************************************************
*  Find index of next bit set in a bitmap
************************************************************************/




static inline
uint32_t find_next_bit_set(uint64_t *bitmap, uint32_t start_index)
// Return index of lowest bit in bitmap given by ``bitmap`` set to one.
// We start searching the bitmap at index ``start_index``.
// Caution: 
// At least one bit in the bitmap with index >= ``start_index`` must
// be set. Otherwise buffer overflow occurs!
{
     uint64_t v, *pb = bitmap + (start_index >> 6);
     uint32_t index;
     v = *pb & (0ULL - ((uint64_t)1ULL << (uint64_t)(start_index & 63)));
     while (v == 0) v = *++pb;
     index = (uint32_t)(pb - bitmap) << 6 ;      
     // Return index + (position of least significant bit of v)
     // Use a deBruijn sequence to locate the least significant bit
     v &= 0ULL - v;
     v = (uint64_t)(v * UINT64T_LOWBIT_MULTIPLIER);
     return index + UINT64T_LOWBIT_TABLE[v >> 57]; 
}

#define get_bitmap(j) (uint32_t)((bitmap[(j) >> 6] >> ((j) & 63)) & 1UL)
#define set_bitmap(j) bitmap[(j) >> 6] |= ((uint64_t)1ULL << ((j) & 63))



/************************************************************************
*  Auxiliary functions manipulating stuctures for orbit arrays
************************************************************************/


static inline uint32_t lin2_generator_size(uint32_t n, uint32_t k)
// Return number of 32-bit integers required to store ``k`` generators
// of a group acting on a vector space over GF(2) of dimension ``n``.
{
    return 2 * (n + 1) * k;
}

static inline int32_t load_lin2_info(uint32_t *a, lin2_type *ps)
// Copy header and pointer to data from an opaque array ``a``
// into the structure of type ``lin2_type`` referred by ``ps``.
// The function returns the status stored in ``ps->status``.
// A negative return value indicates an error.
{
    if (a == NULL) return ERR_GEN_UFIND_INT_LIN2 - 1;
    ps->status = (int32_t)a[0];
    ps->n = a[1];
    ps->n_max_g = a[2];
    ps->n_g = a[3];
    ps->n_orbits = a[4];
    ps->n_vectors = a[5];
    ps->p_t = a + LIN2_LEN_HEADER;
    if (ps->status == LIN2_COMPRESSED) {
        ps->p_o = ps->p_t + ps->n_vectors;
        ps->p_g = ps->p_o + 2 * ps->n_orbits;
   } else {
        ps->p_o = ps->p_t + ((size_t)1UL << ps->n);
        ps->p_g = ps->p_o + ((size_t)1UL << ps->n);
    }
    return ps->status;
}

static inline int32_t store_lin2_info(lin2_type *ps, uint32_t *a)
// Copy header and pointer to data from a structure of
// type ``lin2_type`` referred by ``ps`` to an opaque array ``a``.
// The function returns the status stored in ``ps->status``.
// A negative return value indicates an error.
{
    if (a == NULL) return ERR_GEN_UFIND_INT_LIN2 - 1;
    if ((int32_t)a[0] < 0) return (int32_t)a[0];
    a[0] = (uint32_t)ps->status;
    a[1] = ps->n;
    a[2] = ps->n_max_g;
    a[3] = ps->n_g;
    a[4] = ps->n_orbits;
    a[5] = ps->n_vectors;
    return ps->status;
}

static inline int32_t lin2_error(uint32_t *a, int32_t status)
// If ``status`` is negative and the opaque array ``a`` is not
// in an error status then the function sets the status of ``a``
// to ``status``. The function returns ``status``. 
{
    if (status >= 0) return status;
    if (a != NULL) {
        if ((int32_t)(a[0]) >= 0) a[0] = (uint32_t)(status);
    }
    return status;
}

/*************************************************************************
** Check that a buffer has sufficient length
*************************************************************************/

static inline int32_t
// Check if a pointer ``p_buf`` to a buffer of ``len`` integers of
// type uint32_t is sufficient to store ``min_len`` integers.
// Return 0 if this is the case and a negative value otherwise.
check_out_buf32(uint32_t *p_buf, uint32_t len, uint32_t min_len)
{
     if (p_buf == NULL || min_len & 0x80000000UL || len < min_len)
        return ERR_GEN_UFIND_OUT_SHORT;
     return 0;
}

/************************************************************************
*  Multiply bit vector with bit matrix
************************************************************************/

static inline uint32_t vmatmul(uint32_t v, uint32_t *m)
// Multiply 32-bit vector with 32 times 32 bit matrix v. You may zero
// high bits of v to multiply with matrix with less than 32 rows.
{
    uint32_t w = 0;
    while (v) {
       w ^= *m++ & (0UL - (v & 1UL));
       v >>= 1;
    }
    return w;
}

/************************************************************************
*  Affine operation on GF(2)**n
************************************************************************/

// Let \f$V = GF_2^n\f$. We store an affine mapping \f$V \rightarrow V\f$
// given by \f$v \mapsto v \cdot A + b, v, b \in V\f$, \f$A\f$ an
// \f$n \times n\f$ bit matrix, in an array ``m`` of \f$n + 1\f$
// unsigned 32-bit integers as follows.
// Row ``i``, ``0 <= i < n`` of matrix \f$A\f$ is stored in ``m[i]`` as
// a bit vector; bit vector ``b`` is stored in  ``m[n]``.


// Transform a vector with an affine transformation as given above 
static inline int32_t
vmatmul_aff(uint32_t v, uint32_t *m, uint32_t n)
{
    uint32_t w = m[n], i;
    for (i = 0; i < n; ++i) w ^= *m++ & (0UL - ((v >> i) & 1UL));
    return w & ((1UL << n) - 1);
}



// Invert an affine transformation as given above
static inline int32_t
mat_inverse_aff(uint32_t *m, uint32_t n, uint32_t *m_inv)
// The function computes the inverse of the affine transformation
// on \f$GF_2^n\f$ given by the array ``m``. It stores the inverse
// in the array ``m_inv``.  Here ``m`` and ``m_inv`` are encoded as
// arrays of lengeh ``n + 1`` as described above.
// The function returns 0 if inverse can be computed and -1 if not.
{
    uint64_t a[LIN2_MAX_N + 1];
    uint32_t i, mask = (1UL << n) - 1;
    if (n == 0 || n > LIN2_MAX_N) return ERR_GEN_UFIND_INVERSE;
    for (i = 0; i <= n; ++i) a[i] = m[i] & mask;
    a[n] |= (uint64_t)1ULL << n;
    if (bitmatrix64_inv(a, n + 1)) return ERR_GEN_UFIND_INVERSE;
    for (i = 0; i <= n; ++i) m_inv[i] = (uint32_t)(a[i] & mask);
    return 0;
}



/************************************************************************
* Return pointer to i-th generator in a structure of type lin2_type
************************************************************************/

static inline uint32_t*
// Return pointer to the ``i``-th generator in a structure of
// type ``lin2_type`` referred by pointer ``ps``.
lin2_generator(lin2_type *ps, uint32_t i)
{
    return ps->p_g + (size_t)(i * (ps->n + 1));
}





/************************************************************************
*  Store length information in an array of 24-bit integers
************************************************************************/

/*
Let ``a`` ba a (long) array of 32-bit integers, where only the lower
24 bit of an entry are used. We want to partition this array into
contiguous smaller subarrays without using extra memory. So we may
store length information in the upper 8 bits of an entry.

We may write a length information into the first entries of contiguous
subarray of an array with function ``write_length_info``. We may read
that length information with function  ``read_length_info``.

If a subarray has length 1, we must encode that length information
in the upper 8 bit of the first entry of that array. Note that the
next entry of the array must be kept free to encode the length
information for the next subarray. A larger length information ``n``
may be encoded in the first up to ``n`` entries of the subarray. If
our array has few contiguous subarrays, we want to scan quickly over
(the first entries of) all these subarrays. The functions in this
section encode such a length information for subarrays in a
suitable way.
*/

static inline void
write_length_info(uint32_t *a, uint32_t length)
// Encode the length information ``length`` in the upper 8 bits of
// the entries of the (sub)array ``a``. At most ``length`` entries
// of array ``a`` are required to encode the number ``length``.
// Legal values are ``1 <= length < 0x80000000``.
{
    uint32_t x;
    x = 0x80 | ((length > 0x3f) << 6) | (length & 0x3f);
    *a = (*a & 0x00ffffffUL) | (x << 24);
    ++a;
    length >>= 6;
    while (length) {
        x = ((length > 0x3f) << 6) | (length & 0x3f);
        *a = (*a & 0x00ffffffUL) | (x << 24);
        ++a;
        length >>= 6;
    }
}

static inline int32_t
read_length_info(uint32_t *a)
// Return the length information previously encoded in the
// upper 8 bits of the entries of (sub)array ``a`` with
// function ``write_length_info``.
// The function returns a negative value in case of an error.
{
    uint32_t x, length = 0, sh = 0;
    x = (*a++ >> 24) & 0xff;
    if ((x & 0x80) == 0) return -1;
    length = x & 0x3f;
    for (sh = 6; (sh < 24) && (x & 0x40); sh += 6) {
        x = (*a++ >> 24) & 0xff;
        if (x & 0x80) return ERR_GEN_UFIND_INT_LIN2 - 1;
        length += (x & 0x3f) << sh;
    }
    return length;
}



/************************************************************************
*  Expand a pair of a generator and its inverse for fast processing
************************************************************************/

static inline void
store64_gen(uint64_t *o, uint32_t n, uint32_t *g0, uint32_t *g1)
// Put mask = (1 << n) - 1, and
// g[j] = ((uint64_t)(g1[j] & mask) << 32) + (g1[j] & mask).
// Let G be the bit matrix with row vectors g[0], ..., g[n-1].
// Put byte(k, i) =  (i >> (8 * k)) & 0xff.
// Then we compute o[0x100 * k + i] = (byte[k, i] << 8 * k) (*) G.
// Here (*) is the matrix product of a bit row vector with a bit matrix.
// For n <= 24, output array  o  must have length 0x300.
// In this case we will have:
// v (*) G = o[byte(0, v)] ^ o[0x100 + byte(1, v)] ^ o[0x200 + byte(2, v)].
// We set the entries of  o used for computing v (*) G in this case only.
{
    uint32_t i, jmax, j, j_hi;
    uint64_t mask, g_hi;
    memset(o, 0, 0x300 * sizeof(uint64_t));
    mask = (1UL << n) - 1;
    mask = (mask << 32ULL) + mask;
    o[0] = g1[n];            // enter constant term of affine operation
    o[0] = ((o[0] << 32ULL) + g0[n]) & mask; // the other constant term
    for (i = 0; i < 24; i += 8) {
        if (n > i) {
            g_hi = 0; j_hi = 0;
            jmax = n - i;
            if (jmax > 8) jmax = 8;
            for (j = 1; j < (1UL << jmax); ++j) {
                if ((j & (j - 1)) == 0) { // is j a power of 2?
                    g_hi = *g1++;
                    g_hi = ((g_hi << 32ULL) + *g0++) & mask;
                    j_hi = j;
                }
                o[j] = g_hi ^ o[j - j_hi];
            }
        }
        o += 0x100;
    }
}



/************************************************************************
*  Compute the partition into the orbits
************************************************************************/

/** @brief Compute the partition into the orbits

Let ``S(l_t)`` be the set of integers ``i`` with ``0 <= i < l_t``,
where ``l_t = 1 << dim``. Let a partition of ``S(l_T)`` into orbits
be stored in the array ``table`` (of size ``l_t``) in the same way
as in the description of the array referred by entry ``p_t`` of a
structure of tye ``lin2_type``. Here we refer to the description of
that structure in the header of this file in case of ``status == 2``.

We store the partition of the set ``S(l_t)`` into orbits as a list of
lists in the array ``map``. The structure of the array ``map`` is as
in the description of the array referred by entry ``p_o`` of a
structure of type ``lin2_type``, in case of ``status == 2``.

Parameter ``l_ind`` must be an upper bound for the number of
orbits contained in the array ``table``; this is used for
allocating temporary storage.

The function returns the actual number of orbits stored in the
array ``table`` in case of success, and a negative number in case
of failure.

This function is coded along the lines of
function ``gen_ufind_partition`` in file ``gen_union_find.c``, which
perfoms a similar task.
*/
static inline int32_t
compute_partition(uint32_t *table, uint32_t dim, uint32_t *map, uint32_t l_ind)
{
    uint32_t l_t = 1UL << dim; // length of the table
    uint32_t mask = l_t - 1;
    uint32_t *next, *ind, n, last_n, p, i=0, j, last, fst;
    int32_t status = -1;

    if (dim > 24) return ERR_GEN_UFIND_LIN2_DIM;
    next = malloc(sizeof(uint32_t) * (l_t + 2 * l_ind));
    if (next == NULL) return ERR_GEN_UFIND_MEM;
    ind = next + l_t + l_ind;

    // Here ``next`` is an array of l_t ``l_t + i.``; and ``ind``
    // is an array of length ``i``; with i <= l_ind. Index ``i``
    // is incremented whenever we see a new set of the partition while
    // iterating thtough the ``table``. Eventually, ``i`` will be the
    // number of sets in the partition.

    // We will store the set {x_1, ..., x_m} with index ``j`` in the
    // arrays ``ind`` and ``next`` as shown in the following figure.
    // This peculiar scheme is dictated by the fact that array ``next``
    // is usually much larger than array ``ind``. Note that the original
    // information about the set is copied from ``table`` to ``ind``
    // and ``next``. In the ``table`` the enties ``x_i, i > 1`` contain
    // a pointer to ``x_1``; thus index ``j`` must be reachable from
    // the entry ``x_1`` stored in the array ``next``. Also we want 
    // the first entry ``x_1`` and the last entry ``x_m`` to be
    // reachable from the index ``j`` of the set  {x_1, ..., x_m}.
    
    //
    //  ind:     +--------------------j
    //           |                    ^
    //           |                    |
    //           v                    v
    //  next:   x_m  -->  x_1  -->  j + l_t  -->  x_2 --> x_3
    //           ^                                         |
    //           |____________ ............ _______________|

    for (n = 0; n < l_t; ++n) {
        p = table[n];                         // In case m = 1 we put:
        if ((p >> 24) == 0xfe) {              //   ind:            +----i 
            status = ERR_GEN_UFIND_INT_TABLE-21; //                |    |
            if (i >= l_ind) goto cleanup;     //                   v    |
            ind[i] = next[n] = l_t + i;       //   next: x_1  -->  i + l_t
            next[l_t + i] = n;                //          ^             |
            ++i;                              //          |_____________|
        } else {
            p &= mask;
            // Now ``n`` is the entry to be appended to the set with 
            // index ``j``; and  ``p`` is the first element of that set
            // (corresponding to ``x_1`` in the figure above).
            if ((table[p] >> 24) != 0xfe) {
                status = ERR_GEN_UFIND_INT_TABLE - 23;
                goto cleanup;
            }
            j = next[p] - l_t; // index of the set
            if (j >= i) {
                status = ERR_GEN_UFIND_INT_TABLE - 24;
                goto cleanup;
            }
            last = ind[j];        // old ``x_m`` in the figure above
            if (next[last] != p) {
                status = ERR_GEN_UFIND_INT_TABLE - 25;
                goto cleanup;
            }
            // append ``n`` to the list [x_1, ..., x_m]
            ind[j] = next[last] = n;   
            next[n] = p;          // let ``n`` point to the 1st element       
        }
    }

    // Store the number of the sets in ``l_ind``.
    l_ind = i;

    // Copy the elements of the sets from array ``next`` back to the
    // array ``map``, so the elements of the same set are adjacent.
    // Store length information for these sets in the array ``map``.
    n = 0;
    for (i = 0; i < l_ind; ++i) {
        last_n = n;
        last = ind[i];               // Last element of the set
        map[n++] = fst = next[last]; // store 1st element of the set
        table[fst] = last_n | 0xfe000000UL;  
                                     // Let table entry point to index
        p = next[l_t + i];           // 2nd element of set or end marker
        while (p != fst) {           // Here  ``fst`` will be end marker
            map[n++] = p;            // Fill ``map`` until end marker found
            p = next[p];
        }
        write_length_info(map + last_n, n - last_n);
    }

    status = ERR_GEN_UFIND_INT_TABLE - 26;
    if (n != l_t) goto cleanup;
    status = l_ind;      // return the number of the sets
    
cleanup:
    free(next);
    return status;

} 





/************************************************************************
*  End of header file
************************************************************************/


/// @endcond

#endif // ifndef GEN_UFIND_LIN2_AUX_H


