/////////////////////////////////////////////////////////////////////////////
// This C header file has been created automatically. Do not edit!!!
/////////////////////////////////////////////////////////////////////////////

#ifndef MAT24_XI_FUNCTIONS_H
#define MAT24_XI_FUNCTIONS_H




#ifdef __cplusplus
extern "C" {
#endif
uint32_t mat24_xi_g_gray(uint32_t v1);
uint32_t mat24_xi_w2_gray(uint32_t v1);
uint32_t mat24_xi_g_cocode(uint32_t v1);
uint32_t mat24_xi_w2_cocode(uint32_t v1);
uint32_t mat24_xi_op_leech(uint32_t x1, uint32_t u_exp);
uint32_t mat24_xi_mul_leech(uint32_t x1, uint32_t x2);
uint32_t mat24_xi_leech_to_short(uint32_t x1);
uint32_t mat24_xi_short_to_leech(uint32_t x1);
uint32_t mat24_xi_op_short(uint32_t x1, uint32_t u_exp);
uint32_t mat24_xi_make_table(uint32_t u_box, uint32_t u_exp, uint16_t *ptab);
void mat24_xi_invert_table(uint16_t *ptab, uint32_t len, uint32_t ncols, uint16_t *pres, uint32_t len_res);
void mat24_xi_split_table(uint16_t *ptab, uint32_t len, uint32_t mod, uint32_t *psign);
#ifdef __cplusplus
}
#endif
#endif  // #ifndef MAT24_FUNCTIONS_H




