/** @file gen_ufind_lin2_aux.h

  Internal header file for files ``gen_ufind_lin2*.c``.
  The file is simply copied to iths destnation without 
  preprocessing by the code generator.
*/





#ifndef GEN_UFIND_LIN2_AUX_H
#define GEN_UFIND_LIN2_AUX_H



/*************************************************************************
** Error codes 
*************************************************************************/


/// @cond DO_NOT_DOCUMENT

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
/// @endcond 




/************************************************************************
*  Stuctures and auxiliary functions for orbit arrays
************************************************************************/

/// @cond DO_NOT_DOCUMENT 


// A structure for interpreting the data in the array ``a`` in function
// ``gen_ufind_lin2_init``.
//
// Status is interpreted as follows:
// 0: Initial status, generators of the group may be added.
// 1: No more generators may be added, orbit information is available.
//    E.g. function ``gen_ufind_lin2_n_orbits`` upgrades to status 1 
// 2: A Schreier vector has been computed
//    E.g. function ``gen_ufind_lin2_finalize`` upgrades to status 2
// LIN2_COMPRESSED: Array ``a`` has been compressed
//    E.g. function ``gen_ufind_lin2_compress`` upgrades to this status
// Any negative status indicates an error.
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


static inline int32_t load_lin2_info(uint32_t *a, lin2_type *ps)
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
{
    if (status >= 0) return status;
    if (a != NULL) {
        if ((int32_t)(a[0]) >= 0) a[0] = (uint32_t)(status);
    }
    return status;
}



/// @endcond





#endif // ifndef GEN_UFIND_LIN2_AUX_H


