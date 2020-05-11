/////////////////////////////////////////////////////////////////////////////
// This C file has been created automatically. Do not edit!!!
/////////////////////////////////////////////////////////////////////////////

// %%COMMENT
// TODO: Yet to be documented!!!


#include <string.h>
#include "mat24_functions.h"
#include "mm_op7.h"   



static void invert7_xyz(uint_mmv_t *v_in, uint_mmv_t *v_out)
{
    uint_fast32_t i;
    const uint16_t *p_theta = MAT24_THETA_TABLE;
    
    for (i = 0; i <2048; ++i) {
        uint_mmv_t mask = -((uint_mmv_t)(((p_theta[i] >> 12) & 0x1ULL)));
        mask &= 0x7777777777777777ULL;
        *v_out++ = *v_in++ ^ mask;
        mask &= 0x77777777ULL;
        *v_out++ = *v_in++ ^ mask;
    }
}


// %%EXPORT p
void mm_op7_t(uint_mmv_t *v_in,  uint32_t exp, uint_mmv_t *v_out)
{
    uint_mmv_t i, j, exp1;
 
    exp %= 3;
    if (exp == 0) {
        for (i = 0; i < 15468; ++i) v_out[i] = v_in[i];
        return;
    }
    exp1 = 0x1ULL - (uint_mmv_t)exp;

    // Do off-diagonal part of tags A, B, C
    for (i = 0; i < 48; ++i) {
        // %%MUL_MATRIX_T3 v_in, exp1, v_out

        // This is an automatically generated matrix operation, do not change!
        {
        uint_mmv_t r0, r1, r2, r3;

        // Multiply the vector of integers mod 7 stored in
        // (v_in) by t**e, where t is the 3 times 3 triality
        // matrix [[0, 2,  -2], [1, 1, 1], [1,  -1, -1]] / 2.
        // and e = 1 if exp1 = 0, e = 2 if exp1 = 
        // (uint_mmv_t)(-1). The result is stored in (v_out).
        // 
        // v_in and v_out are pointers of type *uint_mmv_t.
        // Components with tags A, B, C referred by (v_in) 
        // are processed, one integer of type uint_mmv_t
        // for each tag.
        // 
        // 
        // Loading vector from rep 196884x with tags A,B,C
        // to v[0...2]. Here v_in refers to the tag A part. 
        // Negate v[2] if exp1 == -1.
        r0 = (v_in)[0];
        r1 = (v_in)[48];
        r2 = (v_in)[96] ^ ((exp1) & 0x7777777777777777ULL);
        // Vector is now  r(i) for i = 0,1,2
        exp1 = ~(exp1);
        r3 = (r1 + (r2 ^ 0x7777777777777777ULL));
        r1 = (r1 + r2);
        r2 = (r3 & 0x8888888888888888ULL);
        r2 = ((r3 - r2) + (r2 >> 3));
        r3 = (r1 & 0x8888888888888888ULL);
        r1 = ((r1 - r3) + (r3 >> 3));
        r1 = (((r1 & 0x1111111111111111ULL) << 2)
            | ((r1 & 0x6666666666666666ULL) >> 1));
        r2 = (((r2 & 0x1111111111111111ULL) << 2)
            | ((r2 & 0x6666666666666666ULL) >> 1));
        r3 = (r0 + (r2 ^ 0x7777777777777777ULL));
        r0 = (r0 + r2);
        r2 = (r3 & 0x8888888888888888ULL);
        r2 = ((r3 - r2) + (r2 >> 3));
        r3 = (r0 & 0x8888888888888888ULL);
        r0 = ((r0 - r3) + (r3 >> 3));
        // Store vector v[0...2] to rep 196884x with 
        // tags A,B,C. Here v_out refers to the tag A part. 
        // Negate v[2] if exp1 == -1.
        (v_out)[0] = r1;
        (v_out)[48] = r0;
        (v_out)[96]  = r2 ^ ((exp1) & 0x7777777777777777ULL);
        exp1 = ~(exp1);
        // 22 lines of code, 38 operations
        }
        // End of automatically generated matrix operation.
 
        ++v_in; ++v_out;
    }

    v_in -= 48;
    v_out -= 48;
    // Do diagonal part of tags A, B, C
    for (i = 0; i < 24; ++i) {
        // Copy diagonal of A, zero diagonals of B and C
        uint_mmv_t mask = 0x7ULL << ((i << 2) & 63);
        j = (i << 1) + (i >> 4);
        v_out[j] = (v_out[j] & ~mask) | (v_in[j] & mask);
        v_out[j + 48] &= ~mask;
        v_out[j + 96] &= ~mask;
        // Zero slack
        j = ((i + 1) << 1) - 1;
        v_out[j] &= 0x77777777ULL;
        v_out[j + 48] &= 0x77777777ULL;
        v_out[j + 96] &= 0x77777777ULL;  
    }


    // Do tag T
    v_in += MM_OP7_OFS_T;
    v_out +=  MM_OP7_OFS_T;
    for (i = 0; i < 759; ++i) {
        // %%MUL_MATRIX_T64 v_in, exp1, v_out

        // This is an automatically generated matrix operation, do not change!
        {
        uint_mmv_t r0, r1, r2, r3, r4;
        uint_mmv_t r5;

        // Multiply the vector of integers mod 7 stored
        // in (v_in) by t**e, where t is the 64 times 64 
        // triality matrix and e = 1 if exp1 = 0, e = 2 if
        // exp1 = (uint_mmv_t)(-1). The result is stored
        // in (v_out).
        // 
        // Loading vector v from array v_in; multiply v
        // with diagonal matrix if exp1 == -1.
        r0 = v_in[0] ^ ((exp1) & 0x7077707777770ULL);
        r1 = v_in[1] ^ ((exp1) & 0x7000000700070777ULL);
        r2 = v_in[2] ^ ((exp1) & 0x7000000700070777ULL);
        r3 = v_in[3] ^ ((exp1) & 0x7770700070000007ULL);
        // Vector is now  r(i) for i = 0,1,2,3
        exp1 = ~(exp1);
        // Exchange component i with component 63-i if i 
        // has odd parity; fix it if i has even parity.
        r4 = ((r0 & 0x770700770070770ULL)
            | (r1 & 0x7007077007707007ULL));
        r4 = ((r4 << 32) | (r4 >> 32));
        r4 = (((r4 & 0x777700007777ULL) << 16)
            | ((r4 >> 16) & 0x777700007777ULL));
        r4 = (((r4 & 0x77007700770077ULL) << 8)
            | ((r4 >> 8) & 0x77007700770077ULL));
        r4 = (((r4 & 0x707070707070707ULL) << 4)
            | ((r4 >> 4) & 0x707070707070707ULL));
        r5 = ((r2 & 0x7007077007707007ULL)
            | (r3 & 0x770700770070770ULL));
        r5 = ((r5 << 32) | (r5 >> 32));
        r5 = (((r5 & 0x777700007777ULL) << 16)
            | ((r5 >> 16) & 0x777700007777ULL));
        r5 = (((r5 & 0x77007700770077ULL) << 8)
            | ((r5 >> 8) & 0x77007700770077ULL));
        r5 = (((r5 & 0x707070707070707ULL) << 4)
            | ((r5 >> 4) & 0x707070707070707ULL));
        r0 = ((r0 & 0x7007077007707007ULL)
            | (r5 & 0x770700770070770ULL));
        r1 = ((r1 & 0x770700770070770ULL)
            | (r5 & 0x7007077007707007ULL));
        r2 = ((r2 & 0x770700770070770ULL)
            | (r4 & 0x7007077007707007ULL));
        r3 = ((r3 & 0x7007077007707007ULL)
            | (r4 & 0x770700770070770ULL));
        // Butterfly: v[i], v[i+1] = v[i]+v[i+1], v[i]-v[i+1]
        r4 = (((r0 << 4) & 0x7070707070707070ULL)
            | ((r0 & 0x7070707070707070ULL) >> 4));
        r0 = ((r0 ^ 0x7070707070707070ULL) + r4);
        r4 = (r0 & 0x8888888888888888ULL);
        r0 = ((r0 - r4) + (r4 >> 3));
        r4 = (((r1 << 4) & 0x7070707070707070ULL)
            | ((r1 & 0x7070707070707070ULL) >> 4));
        r1 = ((r1 ^ 0x7070707070707070ULL) + r4);
        r4 = (r1 & 0x8888888888888888ULL);
        r1 = ((r1 - r4) + (r4 >> 3));
        r4 = (((r2 << 4) & 0x7070707070707070ULL)
            | ((r2 & 0x7070707070707070ULL) >> 4));
        r2 = ((r2 ^ 0x7070707070707070ULL) + r4);
        r4 = (r2 & 0x8888888888888888ULL);
        r2 = ((r2 - r4) + (r4 >> 3));
        r4 = (((r3 << 4) & 0x7070707070707070ULL)
            | ((r3 & 0x7070707070707070ULL) >> 4));
        r3 = ((r3 ^ 0x7070707070707070ULL) + r4);
        r4 = (r3 & 0x8888888888888888ULL);
        r3 = ((r3 - r4) + (r4 >> 3));
        // Vector is now  r(i) for i = 0,1,2,3
        // Butterfly: v[i], v[i+2] = v[i]+v[i+2], v[i]-v[i+2]
        r4 = (((r0 << 8) & 0x7700770077007700ULL)
            | ((r0 & 0x7700770077007700ULL) >> 8));
        r0 = ((r0 ^ 0x7700770077007700ULL) + r4);
        r4 = (r0 & 0x8888888888888888ULL);
        r0 = ((r0 - r4) + (r4 >> 3));
        r4 = (((r1 << 8) & 0x7700770077007700ULL)
            | ((r1 & 0x7700770077007700ULL) >> 8));
        r1 = ((r1 ^ 0x7700770077007700ULL) + r4);
        r4 = (r1 & 0x8888888888888888ULL);
        r1 = ((r1 - r4) + (r4 >> 3));
        r4 = (((r2 << 8) & 0x7700770077007700ULL)
            | ((r2 & 0x7700770077007700ULL) >> 8));
        r2 = ((r2 ^ 0x7700770077007700ULL) + r4);
        r4 = (r2 & 0x8888888888888888ULL);
        r2 = ((r2 - r4) + (r4 >> 3));
        r4 = (((r3 << 8) & 0x7700770077007700ULL)
            | ((r3 & 0x7700770077007700ULL) >> 8));
        r3 = ((r3 ^ 0x7700770077007700ULL) + r4);
        r4 = (r3 & 0x8888888888888888ULL);
        r3 = ((r3 - r4) + (r4 >> 3));
        // Vector is now  r(i) for i = 0,1,2,3
        // Butterfly: v[i], v[i+4] = v[i]+v[i+4], v[i]-v[i+4]
        r4 = (((r0 << 16) & 0x7777000077770000ULL)
            | ((r0 & 0x7777000077770000ULL) >> 16));
        r0 = ((r0 ^ 0x7777000077770000ULL) + r4);
        r4 = (r0 & 0x8888888888888888ULL);
        r0 = ((r0 - r4) + (r4 >> 3));
        r4 = (((r1 << 16) & 0x7777000077770000ULL)
            | ((r1 & 0x7777000077770000ULL) >> 16));
        r1 = ((r1 ^ 0x7777000077770000ULL) + r4);
        r4 = (r1 & 0x8888888888888888ULL);
        r1 = ((r1 - r4) + (r4 >> 3));
        r4 = (((r2 << 16) & 0x7777000077770000ULL)
            | ((r2 & 0x7777000077770000ULL) >> 16));
        r2 = ((r2 ^ 0x7777000077770000ULL) + r4);
        r4 = (r2 & 0x8888888888888888ULL);
        r2 = ((r2 - r4) + (r4 >> 3));
        r4 = (((r3 << 16) & 0x7777000077770000ULL)
            | ((r3 & 0x7777000077770000ULL) >> 16));
        r3 = ((r3 ^ 0x7777000077770000ULL) + r4);
        r4 = (r3 & 0x8888888888888888ULL);
        r3 = ((r3 - r4) + (r4 >> 3));
        // Vector is now  r(i) for i = 0,1,2,3
        // Butterfly: v[i], v[i+8] = v[i]+v[i+8], v[i]-v[i+8]
        r4 = ((r0 << 32) | (r0 >> 32));
        r0 = ((r0 ^ 0x7777777700000000ULL) + r4);
        r4 = (r0 & 0x8888888888888888ULL);
        r0 = ((r0 - r4) + (r4 >> 3));
        r4 = ((r1 << 32) | (r1 >> 32));
        r1 = ((r1 ^ 0x7777777700000000ULL) + r4);
        r4 = (r1 & 0x8888888888888888ULL);
        r1 = ((r1 - r4) + (r4 >> 3));
        r4 = ((r2 << 32) | (r2 >> 32));
        r2 = ((r2 ^ 0x7777777700000000ULL) + r4);
        r4 = (r2 & 0x8888888888888888ULL);
        r2 = ((r2 - r4) + (r4 >> 3));
        r4 = ((r3 << 32) | (r3 >> 32));
        r3 = ((r3 ^ 0x7777777700000000ULL) + r4);
        r4 = (r3 & 0x8888888888888888ULL);
        r3 = ((r3 - r4) + (r4 >> 3));
        // Vector is now  r(i) for i = 0,1,2,3
        // Butterfly: v[i], v[i+16] = v[i]+v[i+16], v[i]-v[i+16]
        r4 = (r0 + (r1 ^ 0x7777777777777777ULL));
        r0 = (r0 + r1);
        r1 = (r4 & 0x8888888888888888ULL);
        r1 = ((r4 - r1) + (r1 >> 3));
        r4 = (r0 & 0x8888888888888888ULL);
        r0 = ((r0 - r4) + (r4 >> 3));
        r4 = (r2 + (r3 ^ 0x7777777777777777ULL));
        r2 = (r2 + r3);
        r3 = (r4 & 0x8888888888888888ULL);
        r3 = ((r4 - r3) + (r3 >> 3));
        r4 = (r2 & 0x8888888888888888ULL);
        r2 = ((r2 - r4) + (r4 >> 3));
        // Vector is now  r(i) for i = 0,1,2,3
        // Butterfly: v[i], v[i+32] = v[i]+v[i+32], v[i]-v[i+32]
        r4 = (r0 + (r2 ^ 0x7777777777777777ULL));
        r0 = (r0 + r2);
        r2 = (r4 & 0x8888888888888888ULL);
        r2 = ((r4 - r2) + (r2 >> 3));
        r4 = (r0 & 0x8888888888888888ULL);
        r0 = ((r0 - r4) + (r4 >> 3));
        r4 = (r1 + (r3 ^ 0x7777777777777777ULL));
        r1 = (r1 + r3);
        r3 = (r4 & 0x8888888888888888ULL);
        r3 = ((r4 - r3) + (r3 >> 3));
        r4 = (r1 & 0x8888888888888888ULL);
        r1 = ((r1 - r4) + (r4 >> 3));
        // Vector is now  r(i) for i = 0,1,2,3
        // Multiplication by 2**-3 is trivial mod 7
        // Storing vector v to array v_out; multiply v
        // with diagonal matrix if exp1 == -1.
        v_out[0] = r0 ^ ((exp1) & 0x7077707777770ULL);
        v_out[1] = r1 ^ ((exp1) & 0x7000000700070777ULL);
        v_out[2] = r2 ^ ((exp1) & 0x7000000700070777ULL);
        v_out[3] = r3 ^ ((exp1) & 0x7770700070000007ULL);
        exp1 = ~(exp1);
        // 112 lines of code, 284 operations
        }
        // End of automatically generated matrix operation.
 
        v_in += 4;
        v_out += 4;
    }

    // Do tags X, Y, and Z
    {
         uint_mmv_t *pXYin, *pYZin, *pZXin;
         uint_mmv_t *pXYout, *pYZout, *pZXout;
         if (exp1 == 0) {
             pXYin = v_in; 
             pXYout = v_out + 8192;  
             pYZin = v_in + 8192; 
             pYZout = v_out + 4096;  
             pZXin = v_in + 4096; 
             pZXout = v_out; 
         } else {
             pXYout = v_out; 
             pXYin = v_in + 8192;  
             pYZout = v_out + 8192; 
             pYZin = v_in + 4096;  
             pZXout = v_out + 4096; 
             pZXin = v_in; 
         }

         // Map X to Y for t and Y to X for t**2
         for (i = 0; i < 4096; ++i) pXYout[i] = pXYin[i];
         mm7_neg_scalprod_d_i(pXYout);
         
         // Map Y to Z for t and Z to Y for t**2
         invert7_xyz(pYZin, pYZout);
         mm7_neg_scalprod_d_i(pYZout);

         // Map Z to X for t and X to Z for t**2
         invert7_xyz(pZXin, pZXout);
    }
}


