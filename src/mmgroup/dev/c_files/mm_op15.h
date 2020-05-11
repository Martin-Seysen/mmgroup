/////////////////////////////////////////////////////////////////////////////
// This C header file has been created automatically. Do not edit!!!
/////////////////////////////////////////////////////////////////////////////

// This header has been created automatically, do not edit!

#ifndef MM_OP15_H
#define MM_OP15_H

#ifdef __cplusplus
extern "C" {
#endif

#include "mm_basics.h"




// Offsets for tags A,B,C,T,X,Z,Y in the internal representation
// in integers of type uint_mmv_t.
// Compare corresponding offsets in file mm_basics.h
#define MM_OP15_OFS_A  (MM_AUX_OFS_A >> 4)
#define MM_OP15_OFS_B  (MM_AUX_OFS_B >> 4)   
#define MM_OP15_OFS_C  (MM_AUX_OFS_C >> 4)    
#define MM_OP15_OFS_T  (MM_AUX_OFS_T >> 4)  
#define MM_OP15_OFS_X  (MM_AUX_OFS_X >> 4)   
#define MM_OP15_OFS_Z  (MM_AUX_OFS_Z >> 4)   
#define MM_OP15_OFS_Y  (MM_AUX_OFS_Y >> 4)    
#define MM_OP15_OFS_E  (MM_AUX_OFS_E >> 4)    
                                  

void mm15_neg_scalprod_d_i(uint_mmv_t* v);
void mm_op15_do_pi(uint_mmv_t *v_in, mm_sub_op_pi_type *p_op, uint_mmv_t *v_out);
void mm_op15_pi(uint_mmv_t *v_in, uint32_t delta, uint32_t pi, uint_mmv_t *v_out);
void mm_op15_delta(uint_mmv_t *v_in, uint32_t delta, uint_mmv_t *v_out);
uint32_t mm_op15_copy(uint_mmv_t *mv1, uint_mmv_t *mv2);
uint32_t mm_op15_compare(uint_mmv_t *mv1, uint_mmv_t *mv2);
void mm_op15_vector_add(uint_mmv_t *mv1, uint_mmv_t *mv2);
void mm_op15_scalar_mul(int32_t factor, uint_mmv_t *mv1);
void mm_op15_do_xy(uint_mmv_t *v_in, mm_sub_op_xy_type *p_op, uint_mmv_t *v_out);
void mm_op15_xy(uint_mmv_t *v_in, uint32_t f, uint32_t e, uint32_t eps, uint_mmv_t *v_out);
void mm_op15_t(uint_mmv_t *v_in,  uint32_t exp, uint_mmv_t *v_out);
void mm_op15_xi(uint_mmv_t *v_in,  uint32_t exp, uint_mmv_t *v_out);
void mm_op15_group_n(uint_mmv_t *v, uint32_t *g, uint_mmv_t *work);
void mm_op15_word(uint_mmv_t *v, uint32_t *g, int32_t len_g, int32_t e, uint_mmv_t *work);
#ifdef __cplusplus
}
#endif
#endif  // #ifndef MM_OP15_H

