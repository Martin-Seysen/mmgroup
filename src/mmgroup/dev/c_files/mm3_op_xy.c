/////////////////////////////////////////////////////////////////////////////
// This C file has been created automatically. Do not edit!!!
/////////////////////////////////////////////////////////////////////////////


// %%COMMENT
// TODO: Adjust this to new order of basis vectors rep!!!

#include <string.h>
#include "mm_op3.h"   
   




static const uint_mmv_t TABLE_PERM64_LOW[] = {
   // %%TABLE TABLE_PERM64_XY_LOW, uint{INT_BITS}
0x0000000000000000ULL,0x0000000000000000ULL,
0x0ff0f00ff00f0ff0ULL,0xf00f0ff00ff0f00fULL,
0x33cccc33cc3333ccULL,0xcc3333cc33cccc33ULL,
0x3c3c3c3c3c3c3c3cULL,0x3c3c3c3c3c3c3c3cULL,
0x3c3cc3c3c3c33c3cULL,0xc3c33c3c3c3cc3c3ULL,
0x33cc33cc33cc33ccULL,0x33cc33cc33cc33ccULL,
0x0ff00ff00ff00ff0ULL,0x0ff00ff00ff00ff0ULL,
0x0000ffffffff0000ULL,0xffff00000000ffffULL,
0x3cc33cc3c33cc33cULL,0xc33cc33c3cc33cc3ULL,
0x3333cccc3333ccccULL,0x3333cccc3333ccccULL,
0x0f0ff0f00f0ff0f0ULL,0x0f0ff0f00f0ff0f0ULL,
0x00ff00ffff00ff00ULL,0xff00ff0000ff00ffULL,
0x00ffff0000ffff00ULL,0x00ffff0000ffff00ULL,
0x0f0f0f0ff0f0f0f0ULL,0xf0f0f0f00f0f0f0fULL,
0x33333333ccccccccULL,0xcccccccc33333333ULL,
0x3cc3c33c3cc3c33cULL,0x3cc3c33c3cc3c33cULL
}; 


static const uint_mmv_t TABLE_PERM64_HIGH[] = {
   // %%TABLE TABLE_PERM64_XY_HIGH, uint{INT_BITS}
0x0000000000000000ULL,0x0000000000000000ULL,
0x3cc3c33c3cc3c33cULL,0xc33c3cc3c33c3cc3ULL,
0xc33c3cc33cc3c33cULL,0xc33c3cc33cc3c33cULL,
0xffffffff00000000ULL,0x00000000ffffffffULL,
0xffffffffffffffffULL,0xffffffffffffffffULL,
0xc33c3cc3c33c3cc3ULL,0x3cc3c33c3cc3c33cULL,
0x3cc3c33cc33c3cc3ULL,0x3cc3c33cc33c3cc3ULL,
0x00000000ffffffffULL,0xffffffff00000000ULL,
0xc003033f033f3ffcULL,0xfcc0c003c003033fULL,
0xfcc0c0033ffcfcc0ULL,0x3ffcfcc0033f3ffcULL,
0x033f3ffc3ffcfcc0ULL,0x3ffcfcc0fcc0c003ULL,
0x3ffcfcc0033f3ffcULL,0xfcc0c0033ffcfcc0ULL,
0x3ffcfcc0fcc0c003ULL,0x033f3ffc3ffcfcc0ULL,
0x033f3ffcc003033fULL,0xc003033ffcc0c003ULL,
0xfcc0c003c003033fULL,0xc003033f033f3ffcULL,
0xc003033ffcc0c003ULL,0x033f3ffcc003033fULL
}; 


static const uint32_t TABLE24_START[4] = {
   MM_OP3_OFS_X, MM_OP3_OFS_Z, MM_OP3_OFS_Y, MM_OP3_OFS_A
};




// %%EXPORT
void mm_op3_do_xy(uint_mmv_t *v_in, mm_sub_op_xy_type *p_op, uint_mmv_t *v_out)
{
    uint_fast32_t i;

    for (i = 0; i < MM_OP3_OFS_E; ++i) v_out[i] = 0;
    
    // Step 1: do rows with 24 entries 
    {
        uint32_t table24_dest[3];
        // TODO: comment properly!!!!
        for (i = 0; i < 3; ++i) table24_dest[i] = TABLE24_START[i];
        i = (TABLE24_START[1] ^ TABLE24_START[2]) & -(p_op->eps & 1);
        table24_dest[1] ^= i;  table24_dest[2] ^= i; 

        for (i = 0; i < 3; ++i)  {
            uint_mmv_t *p_src = v_in + TABLE24_START[i];
            uint_mmv_t *p_dest = v_out + table24_dest[i];
            uint_fast32_t i1;
            uint_mmv_t a_sign[2][1];
            uint_mmv_t d_xor = p_op->lin_d[i];
            uint8_t *p_sign = p_op->sign_XYZ;
    
            for (i1 = 0; i1 < 1; ++i1) {
                uint_mmv_t x = p_op->lin_i[i] >> (i1 << 5); 
                // %%MMV_UINT_SPREAD x, x
                // Spread bits 0,...,31 of x to the (2-bit long) fields
                // of x. A field of x is set to 0x3 if its 
                // corresponding bit in input x is one and to 0 otherwise.
                x = (x & 0xffffULL) 
                    +  ((x & 0xffff0000ULL) << 16);
                x = (x & 0xff000000ffULL) 
                    +  ((x & 0xff000000ff00ULL) << 8);
                x = (x & 0xf000f000f000fULL) 
                    +  ((x & 0xf000f000f000f0ULL) << 4);
                x = (x & 0x303030303030303ULL) 
                    +  ((x & 0xc0c0c0c0c0c0c0cULL) << 2);
                x = (x & 0x1111111111111111ULL) 
                    +  ((x & 0x2222222222222222ULL) << 1);
                x = (((x) << 2) - (x));
                // Bit spreading done.
                a_sign[0][i1] = x;
                a_sign[1][i1] = x ^ 0xffffffffffffffffULL;
            }
            a_sign[1][0] &= 0xffffffffffffULL;
         
            for (i1 = 0; i1 < 2048; ++i1) {
                uint_mmv_t *ps = p_src + ((i1 ^ d_xor) << 0);
                uint_fast8_t sign = (p_sign[i1] >> i) & 1;
                // %%FOR j in range(V24_INTS_USED)
                p_dest[0] = ps[0] ^ a_sign[sign][0];
                // %%END FOR
                p_dest +=  1;      
            }
        }    
    }    

    // Step 2: do rows with 64 entries, tag T // TODO: comment properly!!!!
    {
        uint_mmv_t *p_src = v_in + MM_OP3_OFS_T;
        uint_mmv_t *p_dest = v_out + MM_OP3_OFS_T;
        uint16_t* p_T =  p_op->s_T;
        for (i = 0; i < 759; ++i) {
            uint_fast16_t ofs_l = *p_T;
            uint_fast16_t ofs_h = (ofs_l & 63) >> 5;
            const uint_mmv_t *ps_h = TABLE_PERM64_HIGH +
                ((ofs_l & 0xf000) >> 11);
            const uint_mmv_t *ps_l = TABLE_PERM64_LOW + 
                ((ofs_l & 0xf00) >> 7);
            ofs_l = (ofs_l << 1) & 0x3fULL;
            // %%FOR j in range(V64_INTS)
            p_dest[0] =  ps_h[0] ^ ps_l[0] ^
               (((p_src[0 ^ ofs_h] >> (0 ^ ofs_l)) & 3) << 0) ^
               (((p_src[0 ^ ofs_h] >> (2 ^ ofs_l)) & 3) << 2) ^
               (((p_src[0 ^ ofs_h] >> (4 ^ ofs_l)) & 3) << 4) ^
               (((p_src[0 ^ ofs_h] >> (6 ^ ofs_l)) & 3) << 6) ^
               (((p_src[0 ^ ofs_h] >> (8 ^ ofs_l)) & 3) << 8) ^
               (((p_src[0 ^ ofs_h] >> (10 ^ ofs_l)) & 3) << 10) ^
               (((p_src[0 ^ ofs_h] >> (12 ^ ofs_l)) & 3) << 12) ^
               (((p_src[0 ^ ofs_h] >> (14 ^ ofs_l)) & 3) << 14) ^
               (((p_src[0 ^ ofs_h] >> (16 ^ ofs_l)) & 3) << 16) ^
               (((p_src[0 ^ ofs_h] >> (18 ^ ofs_l)) & 3) << 18) ^
               (((p_src[0 ^ ofs_h] >> (20 ^ ofs_l)) & 3) << 20) ^
               (((p_src[0 ^ ofs_h] >> (22 ^ ofs_l)) & 3) << 22) ^
               (((p_src[0 ^ ofs_h] >> (24 ^ ofs_l)) & 3) << 24) ^
               (((p_src[0 ^ ofs_h] >> (26 ^ ofs_l)) & 3) << 26) ^
               (((p_src[0 ^ ofs_h] >> (28 ^ ofs_l)) & 3) << 28) ^
               (((p_src[0 ^ ofs_h] >> (30 ^ ofs_l)) & 3) << 30) ^
               (((p_src[0 ^ ofs_h] >> (32 ^ ofs_l)) & 3) << 32) ^
               (((p_src[0 ^ ofs_h] >> (34 ^ ofs_l)) & 3) << 34) ^
               (((p_src[0 ^ ofs_h] >> (36 ^ ofs_l)) & 3) << 36) ^
               (((p_src[0 ^ ofs_h] >> (38 ^ ofs_l)) & 3) << 38) ^
               (((p_src[0 ^ ofs_h] >> (40 ^ ofs_l)) & 3) << 40) ^
               (((p_src[0 ^ ofs_h] >> (42 ^ ofs_l)) & 3) << 42) ^
               (((p_src[0 ^ ofs_h] >> (44 ^ ofs_l)) & 3) << 44) ^
               (((p_src[0 ^ ofs_h] >> (46 ^ ofs_l)) & 3) << 46) ^
               (((p_src[0 ^ ofs_h] >> (48 ^ ofs_l)) & 3) << 48) ^
               (((p_src[0 ^ ofs_h] >> (50 ^ ofs_l)) & 3) << 50) ^
               (((p_src[0 ^ ofs_h] >> (52 ^ ofs_l)) & 3) << 52) ^
               (((p_src[0 ^ ofs_h] >> (54 ^ ofs_l)) & 3) << 54) ^
               (((p_src[0 ^ ofs_h] >> (56 ^ ofs_l)) & 3) << 56) ^
               (((p_src[0 ^ ofs_h] >> (58 ^ ofs_l)) & 3) << 58) ^
               (((p_src[0 ^ ofs_h] >> (60 ^ ofs_l)) & 3) << 60) ^
               (((p_src[0 ^ ofs_h] >> (62 ^ ofs_l)) & 3) << 62);
            p_dest[1] =  ps_h[1] ^ ps_l[1] ^
               (((p_src[1 ^ ofs_h] >> (0 ^ ofs_l)) & 3) << 0) ^
               (((p_src[1 ^ ofs_h] >> (2 ^ ofs_l)) & 3) << 2) ^
               (((p_src[1 ^ ofs_h] >> (4 ^ ofs_l)) & 3) << 4) ^
               (((p_src[1 ^ ofs_h] >> (6 ^ ofs_l)) & 3) << 6) ^
               (((p_src[1 ^ ofs_h] >> (8 ^ ofs_l)) & 3) << 8) ^
               (((p_src[1 ^ ofs_h] >> (10 ^ ofs_l)) & 3) << 10) ^
               (((p_src[1 ^ ofs_h] >> (12 ^ ofs_l)) & 3) << 12) ^
               (((p_src[1 ^ ofs_h] >> (14 ^ ofs_l)) & 3) << 14) ^
               (((p_src[1 ^ ofs_h] >> (16 ^ ofs_l)) & 3) << 16) ^
               (((p_src[1 ^ ofs_h] >> (18 ^ ofs_l)) & 3) << 18) ^
               (((p_src[1 ^ ofs_h] >> (20 ^ ofs_l)) & 3) << 20) ^
               (((p_src[1 ^ ofs_h] >> (22 ^ ofs_l)) & 3) << 22) ^
               (((p_src[1 ^ ofs_h] >> (24 ^ ofs_l)) & 3) << 24) ^
               (((p_src[1 ^ ofs_h] >> (26 ^ ofs_l)) & 3) << 26) ^
               (((p_src[1 ^ ofs_h] >> (28 ^ ofs_l)) & 3) << 28) ^
               (((p_src[1 ^ ofs_h] >> (30 ^ ofs_l)) & 3) << 30) ^
               (((p_src[1 ^ ofs_h] >> (32 ^ ofs_l)) & 3) << 32) ^
               (((p_src[1 ^ ofs_h] >> (34 ^ ofs_l)) & 3) << 34) ^
               (((p_src[1 ^ ofs_h] >> (36 ^ ofs_l)) & 3) << 36) ^
               (((p_src[1 ^ ofs_h] >> (38 ^ ofs_l)) & 3) << 38) ^
               (((p_src[1 ^ ofs_h] >> (40 ^ ofs_l)) & 3) << 40) ^
               (((p_src[1 ^ ofs_h] >> (42 ^ ofs_l)) & 3) << 42) ^
               (((p_src[1 ^ ofs_h] >> (44 ^ ofs_l)) & 3) << 44) ^
               (((p_src[1 ^ ofs_h] >> (46 ^ ofs_l)) & 3) << 46) ^
               (((p_src[1 ^ ofs_h] >> (48 ^ ofs_l)) & 3) << 48) ^
               (((p_src[1 ^ ofs_h] >> (50 ^ ofs_l)) & 3) << 50) ^
               (((p_src[1 ^ ofs_h] >> (52 ^ ofs_l)) & 3) << 52) ^
               (((p_src[1 ^ ofs_h] >> (54 ^ ofs_l)) & 3) << 54) ^
               (((p_src[1 ^ ofs_h] >> (56 ^ ofs_l)) & 3) << 56) ^
               (((p_src[1 ^ ofs_h] >> (58 ^ ofs_l)) & 3) << 58) ^
               (((p_src[1 ^ ofs_h] >> (60 ^ ofs_l)) & 3) << 60) ^
               (((p_src[1 ^ ofs_h] >> (62 ^ ofs_l)) & 3) << 62);
            // %%END FOR
            p_src += 2; 
            p_dest += 2; 
            ++p_T;
        }
    }


    // Step 3: do rows with 24 entries, tags A, B, C // TODO: comment properly!!!!
    {
        uint_mmv_t mask[4];
        uint_mmv_t neg_mask[1];
        uint_mmv_t f = p_op->f_i, ef = p_op->ef_i, eps;
        for (i = 0; i < 1; ++i) {
             mask[i] = f >> (i << 5);
             mask[i + 1] = ef >> (i << 5);
        }
        neg_mask[0] = 0xffffffffffffULL;
        for (i = 0; i < 2; ++i) {
            uint_mmv_t x = mask[i];
            // %%MMV_UINT_SPREAD x, x
            // Spread bits 0,...,31 of x to the (2-bit long) fields
            // of x. A field of x is set to 0x3 if its 
            // corresponding bit in input x is one and to 0 otherwise.
            x = (x & 0xffffULL) 
                +  ((x & 0xffff0000ULL) << 16);
            x = (x & 0xff000000ffULL) 
                +  ((x & 0xff000000ff00ULL) << 8);
            x = (x & 0xf000f000f000fULL) 
                +  ((x & 0xf000f000f000f0ULL) << 4);
            x = (x & 0x303030303030303ULL) 
                +  ((x & 0xc0c0c0c0c0c0c0cULL) << 2);
            x = (x & 0x1111111111111111ULL) 
                +  ((x & 0x2222222222222222ULL) << 1);
            x = (((x) << 2) - (x));
            // Bit spreading done.
            mask[i] = x = x & neg_mask[i & 0];
            mask[i + 2] = x ^ neg_mask[i & 0];
        }

        f <<= 1;
        ef <<= 1;
        eps = -(p_op->eps & 0x1ULL);
        for (i = 0; i < 24; i += 1) {
            uint_mmv_t t, t1, t2;
            // %%FOR j in range(V24_INTS_USED)
            // process uint_mmv_t 0 of row i/1 for tags A, B, C
            t1 = v_in[i + MM_OP3_OFS_A + 0]; 
            t = mask[0 + (f & 0x2ULL)];
            v_out[i + MM_OP3_OFS_A + 0] = t1 ^ t; 
            t1 = v_in[i + MM_OP3_OFS_B + 0]; 
            t2 = v_in[i + MM_OP3_OFS_C + 0];
            t &= (t1 ^ t2);
            t ^= mask[1 + (ef & 0x2ULL)];
            v_out[i + MM_OP3_OFS_B + 0] = t1 ^ t;
            t2 ^= eps & neg_mask[0];
            v_out[i + MM_OP3_OFS_C + 0] = t2 ^ t;
            // %%END FOR
            f >>= 1; ef >>= 1;      
        }
        //yet to be done!!!!
    }



    // If eps is odd: 
    //    negate entries X_d,i with scalar product <d,i> = 1
    if (p_op->eps & 1) mm3_neg_scalprod_d_i(v_out + MM_OP3_OFS_X); 
} 





// %%EXPORT p
void mm_op3_xy(uint_mmv_t *v_in, uint32_t f, uint32_t e, uint32_t eps, uint_mmv_t *v_out)
{
    mm_sub_op_xy_type s_op;
    mm_sub_prep_xy(f, e, eps, &s_op);
    mm_op3_do_xy(v_in, &s_op, v_out);
}


