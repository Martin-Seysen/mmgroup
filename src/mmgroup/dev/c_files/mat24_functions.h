/////////////////////////////////////////////////////////////////////////////
// This C header file has been created automatically. Do not edit!!!
/////////////////////////////////////////////////////////////////////////////

#ifndef MAT24_FUNCTIONS_H
#define MAT24_FUNCTIONS_H

#include <stdint.h>

#define MAT24_DLL  // We want a DLL!!


// Generic helper definitions for shared library support
#if defined _WIN32 || defined __CYGWIN__
  #define MAT24_HELPER_DLL_IMPORT __declspec(dllimport)
  #define MAT24_HELPER_DLL_EXPORT __declspec(dllexport)
  #define MAT24_HELPER_DLL_LOCAL
#else
  #if __GNUC__ >= 4
    #define MAT24_HELPER_DLL_IMPORT __attribute__ ((visibility ("default")))
    #define MAT24_HELPER_DLL_EXPORT __attribute__ ((visibility ("default")))
    #define MAT24_HELPER_DLL_LOCAL  __attribute__ ((visibility ("hidden")))
  #else
    #define MAT24_HELPER_DLL_IMPORT
    #define MAT24_HELPER_DLL_EXPORT
    #define MAT24_HELPER_DLL_LOCAL
  #endif
#endif


// Now we use the generic helper definitions above to define MAT24_API 
// and MAT24_LOCAL.
// MAT24_API is used for the public API symbols. It either DLL imports 
// or DLL exports (or does nothing for static build). 
// MAT24_LOCAL is used for non-api symbols.

#ifdef MAT24_DLL // defined if FOX is compiled as a DLL
  #ifdef MAT24_DLL_EXPORTS // defined if we are building the FOX DLL 
                           // (instead of using it)
    #define MAT24_API MAT24_HELPER_DLL_EXPORT
  #else
    #define MAT24_API MAT24_HELPER_DLL_IMPORT
  #endif // MAT24_DLL_EXPORTS
  #define MAT24_LOCAL MAT24_HELPER_DLL_LOCAL
#else // MAT24_DLL is not defined: this means FOX is a static lib.
  #define MAT24_API
  #define MAT24_LOCAL
#endif // MAT24_DLL



#define MAT24_ORDER 244823040 // Order of Mathieu group Mat24

#ifdef __cplusplus
extern "C" {
#endif
MAT24_API
uint32_t mat24_lsbit24(uint32_t v1);
MAT24_API
uint32_t mat24_bw24(uint32_t v1);
MAT24_API
uint32_t mat24_vect_to_bit_list(uint32_t v1, uint8_t *a_out);
MAT24_API
uint32_t mat24_extract_b24(uint32_t v1, uint32_t u_mask);
MAT24_API
uint32_t mat24_spread_b24(uint32_t v1, uint32_t u_mask);
MAT24_API
extern const uint32_t MAT24_BASIS[24];
MAT24_API
extern const uint32_t MAT24_RECIP_BASIS[24+8];
#define mat24_def_octad_to_gcode(o) (MAT24_OCT_DEC_TABLE[o])
#define mat24_def_gcode_to_octad(v) \
  ((MAT24_OCT_ENC_TABLE[(v) >> 1] >> 1) \
    + 3 * ((v) >> 4) - 11)
MAT24_API
extern const uint16_t MAT24_OCT_DEC_TABLE[759];
MAT24_API
extern const uint8_t MAT24_OCT_ENC_TABLE[2048];
MAT24_API
extern const uint16_t MAT24_THETA_TABLE[];
MAT24_API
uint32_t mat24_vect_to_vintern(uint32_t v1);
MAT24_API
uint32_t mat24_vintern_to_vect(uint32_t v1);
MAT24_API
uint32_t mat24_vect_to_cocode(uint32_t v1);
MAT24_API
uint32_t mat24_gcode_to_vect(uint32_t v1);
MAT24_API
uint32_t mat24_cocode_to_vect(uint32_t c1);
MAT24_API
uint32_t mat24_vect_to_gcode(uint32_t v1);
MAT24_API
uint32_t mat24_gcode_to_octad(uint32_t v1);
MAT24_API
uint32_t mat24_vect_to_octad(uint32_t v1);
MAT24_API
uint32_t mat24_octad_to_gcode(uint32_t u_octad);
MAT24_API
uint32_t mat24_octad_to_vect(uint32_t u_octad);
MAT24_API
uint32_t mat24_cocode_syndrome(uint32_t c1, uint32_t u_tetrad);
MAT24_API
uint32_t mat24_syndrome(uint32_t v1, uint32_t u_tetrad);
MAT24_API
uint32_t mat24_gcode_weight(uint32_t v1);
MAT24_API
uint32_t mat24_cocode_weight(uint32_t c1);
MAT24_API
uint32_t mat24_scalar_prod(uint32_t v1, uint32_t c1);
MAT24_API
uint32_t mat24_suboctad_to_cocode(uint32_t u_sub, uint32_t v1);
MAT24_API
uint32_t mat24_cocode_to_suboctad(uint32_t c1, uint32_t v1);
MAT24_API
uint32_t mat24_suboctad_weight(uint32_t u_sub);
MAT24_API
uint32_t mat24_suboctad_scalar_prod(uint32_t u_sub1, uint32_t u_sub2);
MAT24_API
uint32_t mat24_ploop_theta(uint32_t v1);
MAT24_API
uint32_t mat24_ploop_cocycle(uint32_t v1, uint32_t v2);
MAT24_API
uint32_t mat24_mul_ploop(uint32_t v1, uint32_t v2);
MAT24_API
uint32_t mat24_pow_ploop(uint32_t v1, uint32_t u_exp);
MAT24_API
uint32_t mat24_ploop_comm(uint32_t v1, uint32_t v2);
MAT24_API
uint32_t mat24_ploop_cap(uint32_t v1, uint32_t v2);
MAT24_API
uint32_t mat24_ploop_assoc(uint32_t v1, uint32_t v2, uint32_t v3);
MAT24_API
uint32_t mat24_perm_complete_heptad(uint8_t *p_io);
MAT24_API
uint32_t mat24_perm_check(uint8_t *p1);
MAT24_API
uint32_t mat24_perm_from_heptads(uint8_t *h1, uint8_t *h2, uint8_t *p_out);
MAT24_API
uint32_t mat24_m24num_to_perm(uint32_t u_m24, uint8_t *p_out);
MAT24_API
uint32_t mat24_perm_to_m24num(uint8_t  *p1);
MAT24_API
void mat24_perm_to_matrix(uint8_t  *p1, uint32_t *m_out);
MAT24_API
void mat24_matrix_to_perm(uint32_t *m1, uint8_t *p_out);
MAT24_API
uint32_t mat24_op_vect_perm(uint32_t v1, uint8_t *p1);
MAT24_API
uint32_t mat24_op_gcode_matrix(uint32_t v1, uint32_t *m1);
MAT24_API
uint32_t mat24_op_gcode_perm(uint32_t v1, uint8_t *p1);
MAT24_API
uint32_t mat24_op_cocode_perm(uint32_t c1, uint8_t *p1);
MAT24_API
void mat24_mul_perm(uint8_t *p1, uint8_t *p2, uint8_t *p_out);
MAT24_API
void mat24_inv_perm(uint8_t *p1, uint8_t *p_out);
MAT24_API
void mat24_autpl_set_qform(uint32_t *m_io);
MAT24_API
void mat24_perm_to_autpl(uint32_t c1, uint8_t *p1, uint32_t *m_out);
MAT24_API
void mat24_cocode_to_autpl(uint32_t c1, uint32_t *m_out);
MAT24_API
void mat24_autpl_to_perm(uint32_t *m1, uint8_t  *p_out);
MAT24_API
uint32_t mat24_autpl_to_cocode(uint32_t *m1);
MAT24_API
uint32_t mat24_op_ploop_autpl(uint32_t v1, uint32_t *m1);
MAT24_API
void mat24_mul_autpl(uint32_t *m1, uint32_t *m2, uint32_t *m_out);
MAT24_API
void mat24_inv_autpl(uint32_t *m1, uint32_t *m_out);
MAT24_API
void mat24_perm_to_iautpl(uint32_t c1, uint8_t *p1, uint8_t *p_o, uint32_t *m_o);
MAT24_API
void mat24_perm_to_net(uint8_t *p1, uint32_t *a_out);
MAT24_API
void mat24_op_all_autpl(uint32_t *m1, uint16_t *a_out);
MAT24_API
void mat24_op_all_cocode(uint32_t c1, uint8_t *a_out);
#ifdef __cplusplus
}
#endif
#endif  // #ifndef MAT24_FUNCTIONS_H



