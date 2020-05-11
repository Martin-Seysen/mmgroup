/////////////////////////////////////////////////////////////////////////////
// This C header file has been created automatically. Do not edit!!!
/////////////////////////////////////////////////////////////////////////////

#ifndef MM_BASICS_H
#define MM_BASICS_H

#include <stdint.h>
#include "mat24_functions.h"

// #define MM_BASICS_DLL  // We want a DLL!!


// Generic helper definitions for shared library support
#if defined _WIN32 || defined __CYGWIN__
  #define MM_BASICS_HELPER_DLL_IMPORT __declspec(dllimport)
  #define MM_BASICS_HELPER_DLL_EXPORT __declspec(dllexport)
  #define MM_BASICS_HELPER_DLL_LOCAL
#else
  #if __GNUC__ >= 4
    #define MM_BASICS_HELPER_DLL_IMPORT __attribute__ ((visibility ("default")))
    #define MM_BASICS_HELPER_DLL_EXPORT __attribute__ ((visibility ("default")))
    #define MM_BASICS_HELPER_DLL_LOCAL  __attribute__ ((visibility ("hidden")))
  #else
    #define MM_BASICS_HELPER_DLL_IMPORT
    #define MM_BASICS_HELPER_DLL_EXPORT
    #define MM_BASICS_HELPER_DLL_LOCAL
  #endif
#endif


// Now we use the generic helper definitions above to define MM_BASICS_API 
// and MM_BASICS_LOCAL.
// MM_BASICS_API is used for the public API symbols. It either DLL imports 
// or DLL exports (or does nothing for static build). 
// MM_BASICS_LOCAL is used for non-api symbols.

#ifdef MM_BASICS_DLL // defined if FOX is compiled as a DLL
  #ifdef MM_BASICS_DLL_EXPORTS // defined if we are building the FOX DLL 
                           // (instead of using it)
    #define MM_BASICS_API MM_BASICS_HELPER_DLL_EXPORT
  #else
    #define MM_BASICS_API MM_BASICS_HELPER_DLL_IMPORT
  #endif // MM_BASICS_DLL_EXPORTS
  #define MM_BASICS_LOCAL MM_BASICS_HELPER_DLL_LOCAL
#else // MM_BASICS_DLL is not defined: this means FOX is a static lib.
  #define MM_BASICS_API
  #define MM_BASICS_LOCAL
#endif // MM_BASICS_DLL



typedef uint64_t uint_mmv_t;


#ifdef __cplusplus
extern "C" {
#endif


// Return nonzero value if p is a bad modulus,  
// i.e. not p = 2**k - 1 for some 2 <= k <= 8
#define mm_aux_bad_p(p) (((p) & ((p)+1)) | (((p)-3) & ((0UL-256UL))))

// Offsets for tags A,B,C,T,X,Z,Y in the internal representation
#define MM_AUX_OFS_A       0UL
#define MM_AUX_OFS_B     768UL    //    24*32
#define MM_AUX_OFS_C    1536UL    //  2*24*32
#define MM_AUX_OFS_T    2304UL    //  3*24*32
#define MM_AUX_OFS_X   50880UL    //  MM_AUX_OFS_T +    759*64
#define MM_AUX_OFS_Z  116416UL    //  MM_AUX_OFS_X +   2048*32
#define MM_AUX_OFS_Y  181952UL    //  MM_AUX_OFS_X + 2*2048*32
#define MM_AUX_OFS_E  247488UL    //  MM_AUX_OFS_X + 3*2048*32. i.e
                                  //  total length of internal rep

// Offsets for tags A,B,C,T,X,Z,Y in the external representation
#define MM_AUX_XOFS_A      24UL
#define MM_AUX_XOFS_B     300UL    //  24 + 1*276
#define MM_AUX_XOFS_C     576UL    //  24 + 2*276
#define MM_AUX_XOFS_T     852UL    //  24 + 3*276
#define MM_AUX_XOFS_X   49428UL    //  MM_AUX_XOFS_T +    759*64
#define MM_AUX_XOFS_Z   98580UL    //  MM_AUX_XOFS_X +   2048*24
#define MM_AUX_XOFS_Y  147732UL    //  MM_AUX_XOFS_X + 2*2048*24
#define MM_AUX_XOFS_E  196884UL    //  MM_AUX_XOFS_X + 3*2048*24. i.e
                                   //  total length of external rep


// Tags for labels and values of vectors in the representation space
// A multiple of a unit vector with coordinate 'coord' is encoded
// in the bit fields of a 32-bit integers in the form. 
//   coord (tag, par1, par2) 
#define MM_SPACE_TAG_A      0x2000000
#define MM_SPACE_TAG_B      0x4000000
#define MM_SPACE_TAG_C      0x6000000
#define MM_SPACE_TAG_T      0x8000000
#define MM_SPACE_TAG_X      0xA000000
#define MM_SPACE_TAG_Z      0xC000000
#define MM_SPACE_TAG_Y      0xE000000 
// Mask for all tags:
// Use y = (x & MM_SPACE_MASK_PAR1) << MM_SPACE_SHIFT_PAR1
// to set parameter par1 in y to the value x.
#define MM_SPACE_MASK_TAG     0xE000000 
// Mask and shift factor for parameter par1  
// Use y = (x << MM_SPACE_SHIFT_PAR1) & MM_SPACE_MASK_PAR1
// to set parameter par1 in y to the value x.
#define MM_SPACE_MASK_PAR1    0x1FFC000   
#define MM_SPACE_SHIFT_PAR1          14   
// Mask and shift factor for parameter par12  
// Use y = (x << MM_SPACE_SHIFT_PAR2) & MM_SPACE_MASK_PAR2
// to set parameter par2 in y to the value x.
#define MM_SPACE_MASK_PAR2       0x3F00   
#define MM_SPACE_SHIFT_PAR2           8 
// Mask for coordinate:  
// Use y = x  & MM_SPACE_MASK_COORD
// to set the coordiante in y to the value x.
// Caution: some special routines for modulus p = 2**k - 1
// use only th lowest k bits of the coordinate.
#define MM_SPACE_COORD_PAR1    0x1FFC000   

MM_BASICS_API
void mm_aux_read_mmv1(uint32_t p, uint_mmv_t *mv, uint8_t *b, uint32_t len);
MM_BASICS_API
void mm_aux_read_direct_mmv1(uint32_t p, uint_mmv_t *mv, uint8_t *b, uint32_t len);
MM_BASICS_API
void mm_aux_write_mmv1(uint32_t p, uint8_t *b, uint_mmv_t *mv, uint32_t len);
MM_BASICS_API
void mm_aux_read_mmv24(uint32_t p, uint_mmv_t *mv, uint8_t *b, uint32_t len);
MM_BASICS_API
void mm_aux_write_mmv24(uint32_t p, uint8_t *b, uint_mmv_t *mv, uint32_t len);
MM_BASICS_API
uint8_t mm_aux_get_mmv1(uint32_t p, uint_mmv_t *mv, uint32_t i);
MM_BASICS_API
void mm_aux_put_mmv1(uint32_t p, uint8_t value, uint_mmv_t *mv, uint32_t i);
MM_BASICS_API
uint32_t mm_aux_mmv_size(uint32_t p);
MM_BASICS_API
void mm_aux_zero_mmv(uint32_t p, uint_mmv_t *mv);
MM_BASICS_API
uint8_t mm_aux_get_mmv(uint32_t p, uint_mmv_t *mv, uint32_t i);
MM_BASICS_API
void mm_aux_put_mmv(uint32_t p, uint8_t value, uint_mmv_t *mv, uint32_t i);
MM_BASICS_API
void mm_aux_random_mmv(uint32_t p, uint8_t *seed, uint_mmv_t *mv);
MM_BASICS_API
int32_t mm_aux_reduce_mmv(uint32_t p, uint_mmv_t *mv);
MM_BASICS_API
int32_t mm_aux_check_mmv(uint32_t p, uint_mmv_t *mv);
MM_BASICS_API
void mm_aux_small24_expand(uint8_t *b_src, uint8_t *b_dest);
MM_BASICS_API
void mm_aux_small24_compress(uint8_t *b_src, uint8_t *b_dest);
MM_BASICS_API
void mm_aux_mmv_to_bytes(uint32_t p, uint_mmv_t *mv, uint8_t *b);
MM_BASICS_API
void mm_aux_bytes_to_mmv(uint32_t p, uint8_t *b, uint_mmv_t *mv);
MM_BASICS_API
int32_t mm_aux_mmv_to_sparse(uint32_t p, uint_mmv_t *mv, uint32_t *sp);
MM_BASICS_API
void mm_aux_mmv_extract_sparse(uint32_t p, uint_mmv_t *mv, uint32_t *sp, uint32_t length);
MM_BASICS_API
void mm_aux_mmv_add_sparse(uint32_t p, uint32_t *sp, uint32_t length, uint_mmv_t *mv);
MM_BASICS_API
void mm_aux_mmv_set_sparse(uint32_t p, uint_mmv_t *mv, uint32_t *sp, uint32_t length);
MM_BASICS_API
uint32_t mm_aux_index_extern_to_sparse(uint32_t i);
MM_BASICS_API
void mm_aux_array_extern_to_sparse(uint32_t *a, uint32_t len);
MM_BASICS_API
int32_t mm_aux_index_sparse_to_extern(uint32_t i);
MM_BASICS_API
uint32_t mm_sparse_purge(uint32_t *sp, uint32_t length);
MM_BASICS_API
uint32_t mm_sparse_bitsort(uint32_t *a, uint32_t len, uint32_t mask);
MM_BASICS_API
uint32_t mm_sparse_alloc_reduce(uint32_t p, uint32_t *sp, uint32_t length);
MM_BASICS_API
uint32_t mm_sparse_reduce(uint32_t p, uint32_t *sp, uint32_t length);
MM_BASICS_API
void mm_sparse_mul(uint32_t p, uint32_t *sp, uint32_t length, int32_t factor);
#define MM_RNG_SIZE 266
MM_BASICS_API
void mm_rng_seed(uint8_t *seed, uint32_t no, uint8_t *key, uint32_t len);
MM_BASICS_API
void mm_rng_gen_modp(uint8_t p, uint8_t *seed, uint8_t *out, uint32_t len);
MM_BASICS_API
void mm_rng_gen_bytes(uint8_t *seed, uint8_t *out, uint32_t len);

// Auxiliary structure for mm_sub_op_pi_type
typedef struct {
   uint16_t preimage;
   uint8_t perm[6];
} mm_sub_op_pi64_type;



// Structure used for preparing an operation x_delta * x_pi
// See corresponding comment in file mm_tables.c.
typedef struct {
    uint32_t d;            
    uint32_t pi;
    uint8_t perm[24];
    uint8_t inv_perm[24];
    uint32_t benes_net[9];
    uint16_t tbl_perm24_big[2048+72];
    mm_sub_op_pi64_type tbl_perm64[759];
} mm_sub_op_pi_type;

// Structure used for preparing an operation y_f * x_e *x_eps
// See corresponding comment in file mm_tables.c.
typedef struct {
    uint32_t f;            
    uint32_t e;
    uint32_t eps;
    uint32_t f_i;
    uint32_t ef_i;
    uint32_t lin_i[3];
    uint32_t lin_d[3];
    uint8_t sign_XYZ[2048];
    uint16_t s_T[759];
} mm_sub_op_xy_type;
MM_BASICS_API
extern const uint_mmv_t MM_SUB_INV_2048_TABLE[];
MM_BASICS_API
extern const uint8_t MM_SUB_OCTAD_ELEMENT_TABLE[759*8];
MM_BASICS_API
void  mm_sub_prep_pi(uint32_t delta, uint32_t pi, mm_sub_op_pi_type *p_op);
MM_BASICS_API
void  mm_sub_test_prep_pi_64(uint32_t delta, uint32_t pi, uint32_t *p_tbl);
MM_BASICS_API
void  mm_sub_prep_xy(uint32_t f, uint32_t e, uint32_t eps, mm_sub_op_xy_type *p_op);
MM_BASICS_API
void  mm_sub_test_prep_xy(uint32_t f, uint32_t e, uint32_t eps, uint32_t n, uint32_t *p_tbl);
MM_BASICS_API
void mm_group_n_mul_delta_pi(uint32_t *g, uint32_t delta, uint32_t pi);
MM_BASICS_API
void mm_group_n_mul_inv_delta_pi(uint32_t *g, uint32_t delta, uint32_t pi);
MM_BASICS_API
void mm_group_n_mul_x(uint32_t * g, uint32_t e);
MM_BASICS_API
void mm_group_n_mul_y(uint32_t * g, uint32_t f);
MM_BASICS_API
void mm_group_n_mul_t(uint32_t * g, uint32_t exp);
MM_BASICS_API
void mm_group_n_mul_element(uint32_t *g, uint32_t *g1);
MM_BASICS_API
void mm_group_n_mul_inv_element(uint32_t *g, uint32_t *g1);
MM_BASICS_API
void mm_group_n_clear(uint32_t *g);
MM_BASICS_API
uint32_t mm_group_n_mul_atom(uint32_t *g, uint32_t atom);
MM_BASICS_API
uint32_t mm_group_n_to_word(uint32_t *g, uint32_t *word);
MM_BASICS_API
uint32_t mm_group_split_word_n(uint32_t *word, uint32_t length, uint32_t *g);
MM_BASICS_API
uint32_t mm_group_mul_words(uint32_t *w1, uint32_t l1, uint32_t *w2, uint32_t l2, int32_t e);
// Structure used for referring to tables for operator xi
// See corresponding comment in file mm_tables_xi.c
typedef struct {
   uint16_t *p_perm;
   uint32_t *p_sign;
} mm_sub_table_xi_type;

extern mm_sub_table_xi_type mm_sub_table_xi[2][5];
MM_BASICS_API
uint32_t mm_sub_get_table_xi(uint32_t e, uint32_t i, uint32_t j, uint32_t k);

#ifdef __cplusplus
}
#endif
#endif  // #ifndef MM_BASICS_H
