/////////////////////////////////////////////////////////////////////////////
// This C file has been created automatically. Do not edit!!!
/////////////////////////////////////////////////////////////////////////////

// %%COMMENT
// TODO: Yet to be documented!!!


#include <string.h>
#include "mat24_functions.h"
#include "mm_op127.h"   



static void invert127_xyz(uint_mmv_t *v_in, uint_mmv_t *v_out)
{
    uint_fast32_t i;
    const uint16_t *p_theta = MAT24_THETA_TABLE;
    
    for (i = 0; i <2048; ++i) {
        uint_mmv_t mask = -((uint_mmv_t)(((p_theta[i] >> 12) & 0x1ULL)));
        mask &= 0x7f7f7f7f7f7f7f7fULL;
        *v_out++ = *v_in++ ^ mask;
        *v_out++ = *v_in++ ^ mask;
        *v_out++ = *v_in++ ^ mask;
        *v_out++ = 0;
        ++v_in;
    }
}


// %%EXPORT p
void mm_op127_t(uint_mmv_t *v_in,  uint32_t exp, uint_mmv_t *v_out)
{
    uint_mmv_t i, j, exp1;
 
    exp %= 3;
    if (exp == 0) {
        for (i = 0; i < 30936; ++i) v_out[i] = v_in[i];
        return;
    }
    exp1 = 0x1ULL - (uint_mmv_t)exp;

    // Do off-diagonal part of tags A, B, C
    for (i = 0; i < 24; ++i) {
        for (j = 0; j < 3; ++j) {
            // %%MUL_MATRIX_T3 v_in, exp1, v_out

            // This is an automatically generated matrix operation, do not change!
            {
            uint_mmv_t r0, r1, r2, r3;

            // Multiply the vector of integers mod 127 stored in
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
            r1 = (v_in)[96];
            r2 = (v_in)[192] ^ ((exp1) & 0x7f7f7f7f7f7f7f7fULL);
            // Vector is now  r(i) for i = 0,1,2
            exp1 = ~(exp1);
            r3 = (r1 + (r2 ^ 0x7f7f7f7f7f7f7f7fULL));
            r1 = (r1 + r2);
            r2 = (r3 & 0x8080808080808080ULL);
            r2 = ((r3 - r2) + (r2 >> 7));
            r3 = (r1 & 0x8080808080808080ULL);
            r1 = ((r1 - r3) + (r3 >> 7));
            r1 = (((r1 & 0x101010101010101ULL) << 6)
                | ((r1 & 0x7e7e7e7e7e7e7e7eULL) >> 1));
            r2 = (((r2 & 0x101010101010101ULL) << 6)
                | ((r2 & 0x7e7e7e7e7e7e7e7eULL) >> 1));
            r3 = (r0 + (r2 ^ 0x7f7f7f7f7f7f7f7fULL));
            r0 = (r0 + r2);
            r2 = (r3 & 0x8080808080808080ULL);
            r2 = ((r3 - r2) + (r2 >> 7));
            r3 = (r0 & 0x8080808080808080ULL);
            r0 = ((r0 - r3) + (r3 >> 7));
            // Store vector v[0...2] to rep 196884x with 
            // tags A,B,C. Here v_out refers to the tag A part. 
            // Negate v[2] if exp1 == -1.
            (v_out)[0] = r1;
            (v_out)[96] = r0;
            (v_out)[192]  = r2 ^ ((exp1) & 0x7f7f7f7f7f7f7f7fULL);
            exp1 = ~(exp1);
            // 22 lines of code, 38 operations
            }
            // End of automatically generated matrix operation.
 
            ++v_in; ++v_out;
        }
        ++v_in;
        *v_out++ = 0;
    }

    v_in -= 96;
    v_out -= 96;
    // Do diagonal part of tags A, B, C
    for (i = 0; i < 24; ++i) {
        // Copy diagonal of A, zero diagonals of B and C
        uint_mmv_t mask = 0x7fULL << ((i << 3) & 63);
        j = (i << 2) + (i >> 3);
        v_out[j] = (v_out[j] & ~mask) | (v_in[j] & mask);
        v_out[j + 96] &= ~mask;
        v_out[j + 192] &= ~mask;
    }


    // Do tag T
    v_in += MM_OP127_OFS_T;
    v_out +=  MM_OP127_OFS_T;
    for (i = 0; i < 759; ++i) {
        // %%MUL_MATRIX_T64 v_in, exp1, v_out

        // This is an automatically generated matrix operation, do not change!
        {
        uint_mmv_t r0, r1, r2, r3, r4;
        uint_mmv_t r5, r6, r7, r8;

        // Multiply the vector of integers mod 127 stored
        // in (v_in) by t**e, where t is the 64 times 64 
        // triality matrix and e = 1 if exp1 = 0, e = 2 if
        // exp1 = (uint_mmv_t)(-1). The result is stored
        // in (v_out).
        // 
        // Loading vector v from array v_in; multiply v
        // with diagonal matrix if exp1 == -1.
        r0 = v_in[0] ^ ((exp1) & 0x7f7f7f7f7f7f00ULL);
        r1 = v_in[1] ^ ((exp1) & 0x7f007f7f7fULL);
        r2 = v_in[2] ^ ((exp1) & 0x7f007f7f7fULL);
        r3 = v_in[3] ^ ((exp1) & 0x7f0000000000007fULL);
        r4 = v_in[4] ^ ((exp1) & 0x7f007f7f7fULL);
        r5 = v_in[5] ^ ((exp1) & 0x7f0000000000007fULL);
        r6 = v_in[6] ^ ((exp1) & 0x7f0000000000007fULL);
        r7 = v_in[7] ^ ((exp1) & 0x7f7f7f007f000000ULL);
        // Vector is now  r(i) for i = 0,1,2,3,4,5,6,7
        exp1 = ~(exp1);
        // Exchange component i with component 63-i if i 
        // has odd parity; fix it if i has even parity.
        r8 = ((r0 & 0x7f00007f007f7f00ULL)
            | (r7 & 0x7f7f007f00007fULL));
        r8 = ((r8 << 32) | (r8 >> 32));
        r8 = (((r8 & 0x7f7f00007f7fULL) << 16)
            | ((r8 >> 16) & 0x7f7f00007f7fULL));
        r8 = (((r8 & 0x7f007f007f007fULL) << 8)
            | ((r8 >> 8) & 0x7f007f007f007fULL));
        r0 = ((r0 & 0x7f7f007f00007fULL)
            | (r8 & 0x7f00007f007f7f00ULL));
        r7 = ((r7 & 0x7f00007f007f7f00ULL)
            | (r8 & 0x7f7f007f00007fULL));
        r8 = ((r1 & 0x7f7f007f00007fULL)
            | (r6 & 0x7f00007f007f7f00ULL));
        r8 = ((r8 << 32) | (r8 >> 32));
        r8 = (((r8 & 0x7f7f00007f7fULL) << 16)
            | ((r8 >> 16) & 0x7f7f00007f7fULL));
        r8 = (((r8 & 0x7f007f007f007fULL) << 8)
            | ((r8 >> 8) & 0x7f007f007f007fULL));
        r1 = ((r1 & 0x7f00007f007f7f00ULL)
            | (r8 & 0x7f7f007f00007fULL));
        r6 = ((r6 & 0x7f7f007f00007fULL)
            | (r8 & 0x7f00007f007f7f00ULL));
        r8 = ((r2 & 0x7f7f007f00007fULL)
            | (r5 & 0x7f00007f007f7f00ULL));
        r8 = ((r8 << 32) | (r8 >> 32));
        r8 = (((r8 & 0x7f7f00007f7fULL) << 16)
            | ((r8 >> 16) & 0x7f7f00007f7fULL));
        r8 = (((r8 & 0x7f007f007f007fULL) << 8)
            | ((r8 >> 8) & 0x7f007f007f007fULL));
        r2 = ((r2 & 0x7f00007f007f7f00ULL)
            | (r8 & 0x7f7f007f00007fULL));
        r5 = ((r5 & 0x7f7f007f00007fULL)
            | (r8 & 0x7f00007f007f7f00ULL));
        r8 = ((r3 & 0x7f00007f007f7f00ULL)
            | (r4 & 0x7f7f007f00007fULL));
        r8 = ((r8 << 32) | (r8 >> 32));
        r8 = (((r8 & 0x7f7f00007f7fULL) << 16)
            | ((r8 >> 16) & 0x7f7f00007f7fULL));
        r8 = (((r8 & 0x7f007f007f007fULL) << 8)
            | ((r8 >> 8) & 0x7f007f007f007fULL));
        r3 = ((r3 & 0x7f7f007f00007fULL)
            | (r8 & 0x7f00007f007f7f00ULL));
        r4 = ((r4 & 0x7f00007f007f7f00ULL)
            | (r8 & 0x7f7f007f00007fULL));
        // Butterfly: v[i], v[i+1] = v[i]+v[i+1], v[i]-v[i+1]
        r8 = (((r0 << 8) & 0x7f007f007f007f00ULL)
            | ((r0 & 0x7f007f007f007f00ULL) >> 8));
        r0 = ((r0 ^ 0x7f007f007f007f00ULL) + r8);
        r8 = (r0 & 0x8080808080808080ULL);
        r0 = ((r0 - r8) + (r8 >> 7));
        r8 = (((r1 << 8) & 0x7f007f007f007f00ULL)
            | ((r1 & 0x7f007f007f007f00ULL) >> 8));
        r1 = ((r1 ^ 0x7f007f007f007f00ULL) + r8);
        r8 = (r1 & 0x8080808080808080ULL);
        r1 = ((r1 - r8) + (r8 >> 7));
        r8 = (((r2 << 8) & 0x7f007f007f007f00ULL)
            | ((r2 & 0x7f007f007f007f00ULL) >> 8));
        r2 = ((r2 ^ 0x7f007f007f007f00ULL) + r8);
        r8 = (r2 & 0x8080808080808080ULL);
        r2 = ((r2 - r8) + (r8 >> 7));
        r8 = (((r3 << 8) & 0x7f007f007f007f00ULL)
            | ((r3 & 0x7f007f007f007f00ULL) >> 8));
        r3 = ((r3 ^ 0x7f007f007f007f00ULL) + r8);
        r8 = (r3 & 0x8080808080808080ULL);
        r3 = ((r3 - r8) + (r8 >> 7));
        r8 = (((r4 << 8) & 0x7f007f007f007f00ULL)
            | ((r4 & 0x7f007f007f007f00ULL) >> 8));
        r4 = ((r4 ^ 0x7f007f007f007f00ULL) + r8);
        r8 = (r4 & 0x8080808080808080ULL);
        r4 = ((r4 - r8) + (r8 >> 7));
        r8 = (((r5 << 8) & 0x7f007f007f007f00ULL)
            | ((r5 & 0x7f007f007f007f00ULL) >> 8));
        r5 = ((r5 ^ 0x7f007f007f007f00ULL) + r8);
        r8 = (r5 & 0x8080808080808080ULL);
        r5 = ((r5 - r8) + (r8 >> 7));
        r8 = (((r6 << 8) & 0x7f007f007f007f00ULL)
            | ((r6 & 0x7f007f007f007f00ULL) >> 8));
        r6 = ((r6 ^ 0x7f007f007f007f00ULL) + r8);
        r8 = (r6 & 0x8080808080808080ULL);
        r6 = ((r6 - r8) + (r8 >> 7));
        r8 = (((r7 << 8) & 0x7f007f007f007f00ULL)
            | ((r7 & 0x7f007f007f007f00ULL) >> 8));
        r7 = ((r7 ^ 0x7f007f007f007f00ULL) + r8);
        r8 = (r7 & 0x8080808080808080ULL);
        r7 = ((r7 - r8) + (r8 >> 7));
        // Vector is now  r(i) for i = 0,1,2,3,4,5,6,7
        // Butterfly: v[i], v[i+2] = v[i]+v[i+2], v[i]-v[i+2]
        r8 = (((r0 << 16) & 0x7f7f00007f7f0000ULL)
            | ((r0 & 0x7f7f00007f7f0000ULL) >> 16));
        r0 = ((r0 ^ 0x7f7f00007f7f0000ULL) + r8);
        r8 = (r0 & 0x8080808080808080ULL);
        r0 = ((r0 - r8) + (r8 >> 7));
        r8 = (((r1 << 16) & 0x7f7f00007f7f0000ULL)
            | ((r1 & 0x7f7f00007f7f0000ULL) >> 16));
        r1 = ((r1 ^ 0x7f7f00007f7f0000ULL) + r8);
        r8 = (r1 & 0x8080808080808080ULL);
        r1 = ((r1 - r8) + (r8 >> 7));
        r8 = (((r2 << 16) & 0x7f7f00007f7f0000ULL)
            | ((r2 & 0x7f7f00007f7f0000ULL) >> 16));
        r2 = ((r2 ^ 0x7f7f00007f7f0000ULL) + r8);
        r8 = (r2 & 0x8080808080808080ULL);
        r2 = ((r2 - r8) + (r8 >> 7));
        r8 = (((r3 << 16) & 0x7f7f00007f7f0000ULL)
            | ((r3 & 0x7f7f00007f7f0000ULL) >> 16));
        r3 = ((r3 ^ 0x7f7f00007f7f0000ULL) + r8);
        r8 = (r3 & 0x8080808080808080ULL);
        r3 = ((r3 - r8) + (r8 >> 7));
        r8 = (((r4 << 16) & 0x7f7f00007f7f0000ULL)
            | ((r4 & 0x7f7f00007f7f0000ULL) >> 16));
        r4 = ((r4 ^ 0x7f7f00007f7f0000ULL) + r8);
        r8 = (r4 & 0x8080808080808080ULL);
        r4 = ((r4 - r8) + (r8 >> 7));
        r8 = (((r5 << 16) & 0x7f7f00007f7f0000ULL)
            | ((r5 & 0x7f7f00007f7f0000ULL) >> 16));
        r5 = ((r5 ^ 0x7f7f00007f7f0000ULL) + r8);
        r8 = (r5 & 0x8080808080808080ULL);
        r5 = ((r5 - r8) + (r8 >> 7));
        r8 = (((r6 << 16) & 0x7f7f00007f7f0000ULL)
            | ((r6 & 0x7f7f00007f7f0000ULL) >> 16));
        r6 = ((r6 ^ 0x7f7f00007f7f0000ULL) + r8);
        r8 = (r6 & 0x8080808080808080ULL);
        r6 = ((r6 - r8) + (r8 >> 7));
        r8 = (((r7 << 16) & 0x7f7f00007f7f0000ULL)
            | ((r7 & 0x7f7f00007f7f0000ULL) >> 16));
        r7 = ((r7 ^ 0x7f7f00007f7f0000ULL) + r8);
        r8 = (r7 & 0x8080808080808080ULL);
        r7 = ((r7 - r8) + (r8 >> 7));
        // Vector is now  r(i) for i = 0,1,2,3,4,5,6,7
        // Butterfly: v[i], v[i+4] = v[i]+v[i+4], v[i]-v[i+4]
        r8 = ((r0 << 32) | (r0 >> 32));
        r0 = ((r0 ^ 0x7f7f7f7f00000000ULL) + r8);
        r8 = (r0 & 0x8080808080808080ULL);
        r0 = ((r0 - r8) + (r8 >> 7));
        r8 = ((r1 << 32) | (r1 >> 32));
        r1 = ((r1 ^ 0x7f7f7f7f00000000ULL) + r8);
        r8 = (r1 & 0x8080808080808080ULL);
        r1 = ((r1 - r8) + (r8 >> 7));
        r8 = ((r2 << 32) | (r2 >> 32));
        r2 = ((r2 ^ 0x7f7f7f7f00000000ULL) + r8);
        r8 = (r2 & 0x8080808080808080ULL);
        r2 = ((r2 - r8) + (r8 >> 7));
        r8 = ((r3 << 32) | (r3 >> 32));
        r3 = ((r3 ^ 0x7f7f7f7f00000000ULL) + r8);
        r8 = (r3 & 0x8080808080808080ULL);
        r3 = ((r3 - r8) + (r8 >> 7));
        r8 = ((r4 << 32) | (r4 >> 32));
        r4 = ((r4 ^ 0x7f7f7f7f00000000ULL) + r8);
        r8 = (r4 & 0x8080808080808080ULL);
        r4 = ((r4 - r8) + (r8 >> 7));
        r8 = ((r5 << 32) | (r5 >> 32));
        r5 = ((r5 ^ 0x7f7f7f7f00000000ULL) + r8);
        r8 = (r5 & 0x8080808080808080ULL);
        r5 = ((r5 - r8) + (r8 >> 7));
        r8 = ((r6 << 32) | (r6 >> 32));
        r6 = ((r6 ^ 0x7f7f7f7f00000000ULL) + r8);
        r8 = (r6 & 0x8080808080808080ULL);
        r6 = ((r6 - r8) + (r8 >> 7));
        r8 = ((r7 << 32) | (r7 >> 32));
        r7 = ((r7 ^ 0x7f7f7f7f00000000ULL) + r8);
        r8 = (r7 & 0x8080808080808080ULL);
        r7 = ((r7 - r8) + (r8 >> 7));
        // Vector is now  r(i) for i = 0,1,2,3,4,5,6,7
        // Butterfly: v[i], v[i+8] = v[i]+v[i+8], v[i]-v[i+8]
        r8 = (r0 + (r1 ^ 0x7f7f7f7f7f7f7f7fULL));
        r0 = (r0 + r1);
        r1 = (r8 & 0x8080808080808080ULL);
        r1 = ((r8 - r1) + (r1 >> 7));
        r8 = (r0 & 0x8080808080808080ULL);
        r0 = ((r0 - r8) + (r8 >> 7));
        r8 = (r2 + (r3 ^ 0x7f7f7f7f7f7f7f7fULL));
        r2 = (r2 + r3);
        r3 = (r8 & 0x8080808080808080ULL);
        r3 = ((r8 - r3) + (r3 >> 7));
        r8 = (r2 & 0x8080808080808080ULL);
        r2 = ((r2 - r8) + (r8 >> 7));
        r8 = (r4 + (r5 ^ 0x7f7f7f7f7f7f7f7fULL));
        r4 = (r4 + r5);
        r5 = (r8 & 0x8080808080808080ULL);
        r5 = ((r8 - r5) + (r5 >> 7));
        r8 = (r4 & 0x8080808080808080ULL);
        r4 = ((r4 - r8) + (r8 >> 7));
        r8 = (r6 + (r7 ^ 0x7f7f7f7f7f7f7f7fULL));
        r6 = (r6 + r7);
        r7 = (r8 & 0x8080808080808080ULL);
        r7 = ((r8 - r7) + (r7 >> 7));
        r8 = (r6 & 0x8080808080808080ULL);
        r6 = ((r6 - r8) + (r8 >> 7));
        // Vector is now  r(i) for i = 0,1,2,3,4,5,6,7
        // Butterfly: v[i], v[i+16] = v[i]+v[i+16], v[i]-v[i+16]
        r8 = (r0 + (r2 ^ 0x7f7f7f7f7f7f7f7fULL));
        r0 = (r0 + r2);
        r2 = (r8 & 0x8080808080808080ULL);
        r2 = ((r8 - r2) + (r2 >> 7));
        r8 = (r0 & 0x8080808080808080ULL);
        r0 = ((r0 - r8) + (r8 >> 7));
        r8 = (r1 + (r3 ^ 0x7f7f7f7f7f7f7f7fULL));
        r1 = (r1 + r3);
        r3 = (r8 & 0x8080808080808080ULL);
        r3 = ((r8 - r3) + (r3 >> 7));
        r8 = (r1 & 0x8080808080808080ULL);
        r1 = ((r1 - r8) + (r8 >> 7));
        r8 = (r4 + (r6 ^ 0x7f7f7f7f7f7f7f7fULL));
        r4 = (r4 + r6);
        r6 = (r8 & 0x8080808080808080ULL);
        r6 = ((r8 - r6) + (r6 >> 7));
        r8 = (r4 & 0x8080808080808080ULL);
        r4 = ((r4 - r8) + (r8 >> 7));
        r8 = (r5 + (r7 ^ 0x7f7f7f7f7f7f7f7fULL));
        r5 = (r5 + r7);
        r7 = (r8 & 0x8080808080808080ULL);
        r7 = ((r8 - r7) + (r7 >> 7));
        r8 = (r5 & 0x8080808080808080ULL);
        r5 = ((r5 - r8) + (r8 >> 7));
        // Vector is now  r(i) for i = 0,1,2,3,4,5,6,7
        // Butterfly: v[i], v[i+32] = v[i]+v[i+32], v[i]-v[i+32]
        r8 = (r0 + (r4 ^ 0x7f7f7f7f7f7f7f7fULL));
        r0 = (r0 + r4);
        r4 = (r8 & 0x8080808080808080ULL);
        r4 = ((r8 - r4) + (r4 >> 7));
        r8 = (r0 & 0x8080808080808080ULL);
        r0 = ((r0 - r8) + (r8 >> 7));
        r8 = (r1 + (r5 ^ 0x7f7f7f7f7f7f7f7fULL));
        r1 = (r1 + r5);
        r5 = (r8 & 0x8080808080808080ULL);
        r5 = ((r8 - r5) + (r5 >> 7));
        r8 = (r1 & 0x8080808080808080ULL);
        r1 = ((r1 - r8) + (r8 >> 7));
        r8 = (r2 + (r6 ^ 0x7f7f7f7f7f7f7f7fULL));
        r2 = (r2 + r6);
        r6 = (r8 & 0x8080808080808080ULL);
        r6 = ((r8 - r6) + (r6 >> 7));
        r8 = (r2 & 0x8080808080808080ULL);
        r2 = ((r2 - r8) + (r8 >> 7));
        r8 = (r3 + (r7 ^ 0x7f7f7f7f7f7f7f7fULL));
        r3 = (r3 + r7);
        r7 = (r8 & 0x8080808080808080ULL);
        r7 = ((r8 - r7) + (r7 >> 7));
        r8 = (r3 & 0x8080808080808080ULL);
        r3 = ((r3 - r8) + (r8 >> 7));
        // Vector is now  r(i) for i = 0,1,2,3,4,5,6,7
        // Multiply vector by scalar 2**-3 mod 127
        r0 = (((r0 & 0x707070707070707ULL) << 4)
            | ((r0 & 0x7878787878787878ULL) >> 3));
        r1 = (((r1 & 0x707070707070707ULL) << 4)
            | ((r1 & 0x7878787878787878ULL) >> 3));
        r2 = (((r2 & 0x707070707070707ULL) << 4)
            | ((r2 & 0x7878787878787878ULL) >> 3));
        r3 = (((r3 & 0x707070707070707ULL) << 4)
            | ((r3 & 0x7878787878787878ULL) >> 3));
        r4 = (((r4 & 0x707070707070707ULL) << 4)
            | ((r4 & 0x7878787878787878ULL) >> 3));
        r5 = (((r5 & 0x707070707070707ULL) << 4)
            | ((r5 & 0x7878787878787878ULL) >> 3));
        r6 = (((r6 & 0x707070707070707ULL) << 4)
            | ((r6 & 0x7878787878787878ULL) >> 3));
        r7 = (((r7 & 0x707070707070707ULL) << 4)
            | ((r7 & 0x7878787878787878ULL) >> 3));
        // Storing vector v to array v_out; multiply v
        // with diagonal matrix if exp1 == -1.
        v_out[0] = r0 ^ ((exp1) & 0x7f7f7f7f7f7f00ULL);
        v_out[1] = r1 ^ ((exp1) & 0x7f007f7f7fULL);
        v_out[2] = r2 ^ ((exp1) & 0x7f007f7f7fULL);
        v_out[3] = r3 ^ ((exp1) & 0x7f0000000000007fULL);
        v_out[4] = r4 ^ ((exp1) & 0x7f007f7f7fULL);
        v_out[5] = r5 ^ ((exp1) & 0x7f0000000000007fULL);
        v_out[6] = r6 ^ ((exp1) & 0x7f0000000000007fULL);
        v_out[7] = r7 ^ ((exp1) & 0x7f7f7f007f000000ULL);
        exp1 = ~(exp1);
        // 218 lines of code, 542 operations
        }
        // End of automatically generated matrix operation.
 
        v_in += 8;
        v_out += 8;
    }

    // Do tags X, Y, and Z
    {
         uint_mmv_t *pXYin, *pYZin, *pZXin;
         uint_mmv_t *pXYout, *pYZout, *pZXout;
         if (exp1 == 0) {
             pXYin = v_in; 
             pXYout = v_out + 16384;  
             pYZin = v_in + 16384; 
             pYZout = v_out + 8192;  
             pZXin = v_in + 8192; 
             pZXout = v_out; 
         } else {
             pXYout = v_out; 
             pXYin = v_in + 16384;  
             pYZout = v_out + 16384; 
             pYZin = v_in + 8192;  
             pZXout = v_out + 8192; 
             pZXin = v_in; 
         }

         // Map X to Y for t and Y to X for t**2
         for (i = 0; i < 8192; ++i) pXYout[i] = pXYin[i];
         mm127_neg_scalprod_d_i(pXYout);
         
         // Map Y to Z for t and Z to Y for t**2
         invert127_xyz(pYZin, pYZout);
         mm127_neg_scalprod_d_i(pYZout);

         // Map Z to X for t and X to Z for t**2
         invert127_xyz(pZXin, pZXout);
    }
}


