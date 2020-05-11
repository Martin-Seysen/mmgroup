/////////////////////////////////////////////////////////////////////////////
// This C file has been created automatically. Do not edit!!!
/////////////////////////////////////////////////////////////////////////////

// %%COMMENT
// TODO: Yet to be documented!!!


#include <string.h>
#include "mat24_functions.h"
#include "mm_op3.h"   



static void invert3_xyz(uint_mmv_t *v_in, uint_mmv_t *v_out)
{
    uint_fast32_t i;
    const uint16_t *p_theta = MAT24_THETA_TABLE;
    
    for (i = 0; i <2048; ++i) {
        uint_mmv_t mask = -((uint_mmv_t)(((p_theta[i] >> 12) & 0x1ULL)));
        mask &= 0xffffffffffffULL;
        *v_out++ = *v_in++ ^ mask;
    }
}


// %%EXPORT p
void mm_op3_t(uint_mmv_t *v_in,  uint32_t exp, uint_mmv_t *v_out)
{
    uint_mmv_t i, j, exp1;
 
    exp %= 3;
    if (exp == 0) {
        for (i = 0; i < 7734; ++i) v_out[i] = v_in[i];
        return;
    }
    exp1 = 0x1ULL - (uint_mmv_t)exp;

    // Do off-diagonal part of tags A, B, C
    for (i = 0; i < 24; ++i) {
        // %%MUL_MATRIX_T3 v_in, exp1, v_out

        // This is an automatically generated matrix operation, do not change!
        {
        uint_mmv_t r0, r1, r2, r3, r4;
        uint_mmv_t r5, r6;

        // Multiply the vector of integers mod 3 stored in
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
        r1 = (v_in)[24];
        r2 = (v_in)[48] ^ ((exp1) & 0xffffffffffffffffULL);
        // Vector is now  r(i) for i = 0,1,2,3,4,5
        exp1 = ~(exp1);
        r3 = ((r0 >> 2) & 0x3333333333333333ULL);
        r0 = (r0 & 0x3333333333333333ULL);
        r4 = ((r1 >> 2) & 0x3333333333333333ULL);
        r1 = (r1 & 0x3333333333333333ULL);
        r5 = ((r2 >> 2) & 0x3333333333333333ULL);
        r2 = (r2 & 0x3333333333333333ULL);
        r6 = (r4 + (r5 ^ 0x3333333333333333ULL));
        r4 = (r4 + r5);
        r5 = (r6 & 0x4444444444444444ULL);
        r5 = ((r6 - r5) + (r5 >> 2));
        r6 = (r4 & 0x4444444444444444ULL);
        r4 = ((r4 - r6) + (r6 >> 2));
        r4 = (((r4 & 0x5555555555555555ULL) << 1)
            | ((r4 & 0xaaaaaaaaaaaaaaaaULL) >> 1));
        r5 = (((r5 & 0x5555555555555555ULL) << 1)
            | ((r5 & 0xaaaaaaaaaaaaaaaaULL) >> 1));
        r6 = (r3 + (r5 ^ 0x3333333333333333ULL));
        r3 = (r3 + r5);
        r5 = (r6 & 0x4444444444444444ULL);
        r5 = ((r6 - r5) + (r5 >> 2));
        r6 = (r3 & 0x4444444444444444ULL);
        r3 = ((r3 - r6) + (r6 >> 2));
        r6 = (r1 + (r2 ^ 0x3333333333333333ULL));
        r1 = (r1 + r2);
        r2 = (r6 & 0x4444444444444444ULL);
        r2 = ((r6 - r2) + (r2 >> 2));
        r6 = (r1 & 0x4444444444444444ULL);
        r1 = ((r1 - r6) + (r6 >> 2));
        r1 = (((r1 & 0x5555555555555555ULL) << 1)
            | ((r1 & 0xaaaaaaaaaaaaaaaaULL) >> 1));
        r2 = (((r2 & 0x5555555555555555ULL) << 1)
            | ((r2 & 0xaaaaaaaaaaaaaaaaULL) >> 1));
        r6 = (r0 + (r2 ^ 0x3333333333333333ULL));
        r0 = (r0 + r2);
        r2 = (r6 & 0x4444444444444444ULL);
        r2 = ((r6 - r2) + (r2 >> 2));
        r6 = (r0 & 0x4444444444444444ULL);
        r0 = ((r0 - r6) + (r6 >> 2));
        r0 ^= (r3 << 2);
        r1 ^= (r4 << 2);
        r2 ^= (r5 << 2);
        // Store vector v[0...2] to rep 196884x with 
        // tags A,B,C. Here v_out refers to the tag A part. 
        // Negate v[2] if exp1 == -1.
        (v_out)[0] = r1;
        (v_out)[24] = r0;
        (v_out)[48]  = r2 ^ ((exp1) & 0xffffffffffffffffULL);
        exp1 = ~(exp1);
        // 45 lines of code, 85 operations
        }
        // End of automatically generated matrix operation.
 
        ++v_in; ++v_out;
    }

    v_in -= 24;
    v_out -= 24;
    // Do diagonal part of tags A, B, C
    for (i = 0; i < 24; ++i) {
        // Copy diagonal of A, zero diagonals of B and C
        uint_mmv_t mask = 0x3ULL << ((i << 1) & 63);
        j = (i << 0) + (i >> 5);
        v_out[j] = (v_out[j] & ~mask) | (v_in[j] & mask);
        v_out[j + 24] &= ~mask;
        v_out[j + 48] &= ~mask;
        // Zero slack
        j = ((i + 1) << 0) - 1;
        v_out[j] &= 0xffffffffffffULL;
        v_out[j + 24] &= 0xffffffffffffULL;
        v_out[j + 48] &= 0xffffffffffffULL;  
    }


    // Do tag T
    v_in += MM_OP3_OFS_T;
    v_out +=  MM_OP3_OFS_T;
    for (i = 0; i < 759; ++i) {
        // %%MUL_MATRIX_T64 v_in, exp1, v_out

        // This is an automatically generated matrix operation, do not change!
        {
        uint_mmv_t r0, r1, r2, r3, r4;

        // Multiply the vector of integers mod 3 stored
        // in (v_in) by t**e, where t is the 64 times 64 
        // triality matrix and e = 1 if exp1 = 0, e = 2 if
        // exp1 = (uint_mmv_t)(-1). The result is stored
        // in (v_out).
        // 
        // Loading vector v from array v_in; multiply v
        // with diagonal matrix if exp1 == -1.
        r0 = v_in[0] ^ ((exp1) & 0xc003033f033f3ffcULL);
        r1 = v_in[1] ^ ((exp1) & 0xfcc0c003c003033fULL);
        // Vector is now  r(i) for i = 0,1
        exp1 = ~(exp1);
        // Exchange component i with component 63-i if i 
        // has odd parity; fix it if i has even parity.
        r2 = ((r0 & 0xc33c3cc33cc3c33cULL)
            | (r1 & 0x3cc3c33cc33c3cc3ULL));
        r2 = ((r2 << 32) | (r2 >> 32));
        r2 = (((r2 & 0xffff0000ffffULL) << 16)
            | ((r2 >> 16) & 0xffff0000ffffULL));
        r2 = (((r2 & 0xff00ff00ff00ffULL) << 8)
            | ((r2 >> 8) & 0xff00ff00ff00ffULL));
        r2 = (((r2 & 0xf0f0f0f0f0f0f0fULL) << 4)
            | ((r2 >> 4) & 0xf0f0f0f0f0f0f0fULL));
        r2 = (((r2 & 0x3333333333333333ULL) << 2)
            | ((r2 >> 2) & 0x3333333333333333ULL));
        r0 = ((r0 & 0x3cc3c33cc33c3cc3ULL)
            | (r2 & 0xc33c3cc33cc3c33cULL));
        r1 = ((r1 & 0xc33c3cc33cc3c33cULL)
            | (r2 & 0x3cc3c33cc33c3cc3ULL));
        // Expansion for Hadamard operation:
        // There is no space for a carry bit between bit fields. So 
        // we move bit field 2*i + 1  to bit field 2*i + 64.
        r2 = ((r0 >> 2) & 0x3333333333333333ULL);
        r0 = (r0 & 0x3333333333333333ULL);
        r3 = ((r1 >> 2) & 0x3333333333333333ULL);
        r1 = (r1 & 0x3333333333333333ULL);
        // Vector is now  r(i) for i = 0,1,2,3
        // Butterfly: v[i], v[i+2] = v[i]+v[i+2], v[i]-v[i+2]
        r4 = (((r0 << 4) & 0x3030303030303030ULL)
            | ((r0 & 0x3030303030303030ULL) >> 4));
        r0 = ((r0 ^ 0x3030303030303030ULL) + r4);
        r4 = (r0 & 0x4444444444444444ULL);
        r0 = ((r0 - r4) + (r4 >> 2));
        r4 = (((r1 << 4) & 0x3030303030303030ULL)
            | ((r1 & 0x3030303030303030ULL) >> 4));
        r1 = ((r1 ^ 0x3030303030303030ULL) + r4);
        r4 = (r1 & 0x4444444444444444ULL);
        r1 = ((r1 - r4) + (r4 >> 2));
        r4 = (((r2 << 4) & 0x3030303030303030ULL)
            | ((r2 & 0x3030303030303030ULL) >> 4));
        r2 = ((r2 ^ 0x3030303030303030ULL) + r4);
        r4 = (r2 & 0x4444444444444444ULL);
        r2 = ((r2 - r4) + (r4 >> 2));
        r4 = (((r3 << 4) & 0x3030303030303030ULL)
            | ((r3 & 0x3030303030303030ULL) >> 4));
        r3 = ((r3 ^ 0x3030303030303030ULL) + r4);
        r4 = (r3 & 0x4444444444444444ULL);
        r3 = ((r3 - r4) + (r4 >> 2));
        // Vector is now  r(i) for i = 0,1,2,3
        // Butterfly: v[i], v[i+4] = v[i]+v[i+4], v[i]-v[i+4]
        r4 = (((r0 << 8) & 0x3300330033003300ULL)
            | ((r0 & 0x3300330033003300ULL) >> 8));
        r0 = ((r0 ^ 0x3300330033003300ULL) + r4);
        r4 = (r0 & 0x4444444444444444ULL);
        r0 = ((r0 - r4) + (r4 >> 2));
        r4 = (((r1 << 8) & 0x3300330033003300ULL)
            | ((r1 & 0x3300330033003300ULL) >> 8));
        r1 = ((r1 ^ 0x3300330033003300ULL) + r4);
        r4 = (r1 & 0x4444444444444444ULL);
        r1 = ((r1 - r4) + (r4 >> 2));
        r4 = (((r2 << 8) & 0x3300330033003300ULL)
            | ((r2 & 0x3300330033003300ULL) >> 8));
        r2 = ((r2 ^ 0x3300330033003300ULL) + r4);
        r4 = (r2 & 0x4444444444444444ULL);
        r2 = ((r2 - r4) + (r4 >> 2));
        r4 = (((r3 << 8) & 0x3300330033003300ULL)
            | ((r3 & 0x3300330033003300ULL) >> 8));
        r3 = ((r3 ^ 0x3300330033003300ULL) + r4);
        r4 = (r3 & 0x4444444444444444ULL);
        r3 = ((r3 - r4) + (r4 >> 2));
        // Vector is now  r(i) for i = 0,1,2,3
        // Butterfly: v[i], v[i+8] = v[i]+v[i+8], v[i]-v[i+8]
        r4 = (((r0 << 16) & 0x3333000033330000ULL)
            | ((r0 & 0x3333000033330000ULL) >> 16));
        r0 = ((r0 ^ 0x3333000033330000ULL) + r4);
        r4 = (r0 & 0x4444444444444444ULL);
        r0 = ((r0 - r4) + (r4 >> 2));
        r4 = (((r1 << 16) & 0x3333000033330000ULL)
            | ((r1 & 0x3333000033330000ULL) >> 16));
        r1 = ((r1 ^ 0x3333000033330000ULL) + r4);
        r4 = (r1 & 0x4444444444444444ULL);
        r1 = ((r1 - r4) + (r4 >> 2));
        r4 = (((r2 << 16) & 0x3333000033330000ULL)
            | ((r2 & 0x3333000033330000ULL) >> 16));
        r2 = ((r2 ^ 0x3333000033330000ULL) + r4);
        r4 = (r2 & 0x4444444444444444ULL);
        r2 = ((r2 - r4) + (r4 >> 2));
        r4 = (((r3 << 16) & 0x3333000033330000ULL)
            | ((r3 & 0x3333000033330000ULL) >> 16));
        r3 = ((r3 ^ 0x3333000033330000ULL) + r4);
        r4 = (r3 & 0x4444444444444444ULL);
        r3 = ((r3 - r4) + (r4 >> 2));
        // Vector is now  r(i) for i = 0,1,2,3
        // Butterfly: v[i], v[i+16] = v[i]+v[i+16], v[i]-v[i+16]
        r4 = ((r0 << 32) | (r0 >> 32));
        r0 = ((r0 ^ 0x3333333300000000ULL) + r4);
        r4 = (r0 & 0x4444444444444444ULL);
        r0 = ((r0 - r4) + (r4 >> 2));
        r4 = ((r1 << 32) | (r1 >> 32));
        r1 = ((r1 ^ 0x3333333300000000ULL) + r4);
        r4 = (r1 & 0x4444444444444444ULL);
        r1 = ((r1 - r4) + (r4 >> 2));
        r4 = ((r2 << 32) | (r2 >> 32));
        r2 = ((r2 ^ 0x3333333300000000ULL) + r4);
        r4 = (r2 & 0x4444444444444444ULL);
        r2 = ((r2 - r4) + (r4 >> 2));
        r4 = ((r3 << 32) | (r3 >> 32));
        r3 = ((r3 ^ 0x3333333300000000ULL) + r4);
        r4 = (r3 & 0x4444444444444444ULL);
        r3 = ((r3 - r4) + (r4 >> 2));
        // Vector is now  r(i) for i = 0,1,2,3
        // Butterfly: v[i], v[i+32] = v[i]+v[i+32], v[i]-v[i+32]
        r4 = (r0 + (r1 ^ 0x3333333333333333ULL));
        r0 = (r0 + r1);
        r1 = (r4 & 0x4444444444444444ULL);
        r1 = ((r4 - r1) + (r1 >> 2));
        r4 = (r0 & 0x4444444444444444ULL);
        r0 = ((r0 - r4) + (r4 >> 2));
        r4 = (r2 + (r3 ^ 0x3333333333333333ULL));
        r2 = (r2 + r3);
        r3 = (r4 & 0x4444444444444444ULL);
        r3 = ((r4 - r3) + (r3 >> 2));
        r4 = (r2 & 0x4444444444444444ULL);
        r2 = ((r2 - r4) + (r4 >> 2));
        // Vector is now  r(i) for i = 0,1,2,3
        // Butterfly: v[i], v[i+64] = v[i]+v[i+64], v[i]-v[i+64]
        r4 = (r0 + (r2 ^ 0x3333333333333333ULL));
        r0 = (r0 + r2);
        r2 = (r4 & 0x4444444444444444ULL);
        r2 = ((r4 - r2) + (r2 >> 2));
        r4 = (r0 & 0x4444444444444444ULL);
        r0 = ((r0 - r4) + (r4 >> 2));
        r4 = (r1 + (r3 ^ 0x3333333333333333ULL));
        r1 = (r1 + r3);
        r3 = (r4 & 0x4444444444444444ULL);
        r3 = ((r4 - r3) + (r3 >> 2));
        r4 = (r1 & 0x4444444444444444ULL);
        r1 = ((r1 - r4) + (r4 >> 2));
        // Vector is now  r(i) for i = 0,1,2,3
        // Reverse expansion for Hadamard operation
        r0 ^= (r2 << 2);
        r1 ^= (r3 << 2);
        // Vector is now  r(i) for i = 0,1
        // Multiply vector by scalar 2**-3 mod 3
        r0 = (((r0 & 0x5555555555555555ULL) << 1)
            | ((r0 & 0xaaaaaaaaaaaaaaaaULL) >> 1));
        r1 = (((r1 & 0x5555555555555555ULL) << 1)
            | ((r1 & 0xaaaaaaaaaaaaaaaaULL) >> 1));
        // Storing vector v to array v_out; multiply v
        // with diagonal matrix if exp1 == -1.
        v_out[0] = r0 ^ ((exp1) & 0xc003033f033f3ffcULL);
        v_out[1] = r1 ^ ((exp1) & 0xfcc0c003c003033fULL);
        exp1 = ~(exp1);
        // 110 lines of code, 274 operations
        }
        // End of automatically generated matrix operation.
 
        v_in += 2;
        v_out += 2;
    }

    // Do tags X, Y, and Z
    {
         uint_mmv_t *pXYin, *pYZin, *pZXin;
         uint_mmv_t *pXYout, *pYZout, *pZXout;
         if (exp1 == 0) {
             pXYin = v_in; 
             pXYout = v_out + 4096;  
             pYZin = v_in + 4096; 
             pYZout = v_out + 2048;  
             pZXin = v_in + 2048; 
             pZXout = v_out; 
         } else {
             pXYout = v_out; 
             pXYin = v_in + 4096;  
             pYZout = v_out + 4096; 
             pYZin = v_in + 2048;  
             pZXout = v_out + 2048; 
             pZXin = v_in; 
         }

         // Map X to Y for t and Y to X for t**2
         for (i = 0; i < 2048; ++i) pXYout[i] = pXYin[i];
         mm3_neg_scalprod_d_i(pXYout);
         
         // Map Y to Z for t and Z to Y for t**2
         invert3_xyz(pYZin, pYZout);
         mm3_neg_scalprod_d_i(pYZout);

         // Map Z to X for t and X to Z for t**2
         invert3_xyz(pZXin, pZXout);
    }
}


