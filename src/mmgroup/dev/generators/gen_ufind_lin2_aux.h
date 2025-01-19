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
and stores the relevant information in an array ``a``, which is usually
the first parameter of such a function. This makes the interface to
the outside world simple; but it slightly complicates the internal
operation of a function in module ``gen_ufind_lin2.c``.

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
  1: No more generators may be added, orbit information is available.
     E.g. function ``gen_ufind_lin2_n_orbits`` upgrades to status 1
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

  In case ``0 <= status <= 2`` this table has length ``1 << s.n``,
  where ``s.n`` is the dimension of the vector space \f$V\f$.
  Entry ``i`` in that table corresponds to vector ``i`` in \f$V\f$.
  The union-find algorithm in module ``gen_union_find.c`` is executed
  on this table in order to unite each vector in \f$V\f$ with its
  images under the operation of each generator of the group \f$G\f$.

  Bits 0,...,23 of entry ``i`` of array ``t`` is an index ``j``
  referring to table ``map`` referred by ``s.p_o``. Entry ``map[j]``
  describes the orbit of entry ``i``, as explained in the description
  of entry ``p_o``. This is valid if the status is 1 or 2.

  Bits 24,...,31 of entry ``i`` of array ``t`` is the number ``k`` of
  a generator ``g`` of the group \f$G\f$. Here a Schreier vector is a
  mapping of the vectors in  \f$V\f$ to the union of set of the
  generators of the group \f$G\f$ with the neutral element of \f$G\f$.
  The neutral element of \f$G\f$ is encoded as ``k = 0xfe``.
  This is valid if the status 2.

Structure ``s`` contains an entry ``p_o`` pointing to a table ``map``.
 
  In case ``1 <= status <= 2`` this table has length ``1 << s.n``. It
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


If ``s.status`` is LIN2_COMPRESSED then entries ``s.p_t`` and ``s.p_o``
contain different data. Details will be explained in a future version!


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
        ps->p_o = ps->p_t + 2 * ps->n_orbits;
        ps->p_g = ps->p_o + ps->n_vectors;
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
*  Invert a bit matrix
************************************************************************/

static inline int32_t
mat_inverse(uint32_t *m, uint32_t n, uint32_t *m_inv)
// Store inverse of the ``n`` times ``n`` bit matrix ``m`` in the
// array ``m_inv``. Return 0 if inverse can be computed and -1 if not.
{
    uint64_t a[LIN2_MAX_N];
    uint32_t i, mask = (1UL << n) - 1;
    if (n == 0 || n > LIN2_MAX_N) return ERR_GEN_UFIND_INVERSE;
    for (i = 0; i < n; ++i) a[i] = m[i];
    if (bitmatrix64_inv(a, n)) return ERR_GEN_UFIND_INVERSE;
    for (i = 0; i < n; ++i) m_inv[i] = (uint32_t)(a[i] & mask);
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
* Unite vector v with vector A * v, for a matrix A over GF(2)
************************************************************************/

#define MAT_BLOCKSIZE 7



/** @brief Perform union-find algorithm in ``GF(2)^n``

Let ``S(n)`` be the set of integers ``v`` with ``0 <= v < 1 << n``;
and let a partition of ``S(n)`` be stored in the array ``table``
(of size ``1 << n``) as described in function ``gen_ufind_init``.

In this function the entries of ``S(n)`` are interpreted as bit
vectors. The function joins the set containing ``v`` with the set
containing ``v * g`` for all ``v`` in ``S(n)``. Here ``g`` is an
``n`` times ``n`` bit matrix over GF(2) stored in the array referred
by ``g``. Row ``j`` of bit matrix ``g`` is stored in ``g[j]`` as
an integer encoding a bit vector. All bit vector arithmetic is done
over GF(2).

Thus the array referred by ``table`` must have length ``1 << n``;
and the array referred by ``g`` must have length ``len_g * n``.

The function returns 0 in case of success and -1 in case of error.
*/
static inline int32_t
union_linear(uint32_t *table, uint32_t n, uint32_t *g)
{
     uint32_t j0, j1, lg_bl, bl, w, a[1UL << MAT_BLOCKSIZE];
     uint32_t t_length = 1UL << n;
     uint32_t mask = t_length - 1;
     uint32_t status = 0;

     if (n > LIN2_MAX_N || n == 0) return ERR_GEN_UFIND_LIN2_GEN;
     lg_bl = (n + 1) >> 1;
     lg_bl = lg_bl < MAT_BLOCKSIZE ? lg_bl : MAT_BLOCKSIZE;
     bl = 1UL << lg_bl;
     for (j1 = 0; j1 < bl; ++j1) a[j1] = vmatmul(j1, g) & mask;
     for (j0 = 0; j0 < t_length; j0 += bl) {
         w = vmatmul(j0 >> lg_bl, g + lg_bl) & mask;
         for (j1 = 0; j1 < bl; ++j1) {
             status |=
                gen_ufind_union(table, t_length, j0 ^ j1, w ^ a[j1]);
         }
     }
     return (status & 0x80000000UL) ? ERR_GEN_UFIND_INT_LIN2 - 7 : status;
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



/// @endcond

#endif // ifndef GEN_UFIND_LIN2_AUX_H


