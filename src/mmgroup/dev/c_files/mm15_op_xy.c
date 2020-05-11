/////////////////////////////////////////////////////////////////////////////
// This C file has been created automatically. Do not edit!!!
/////////////////////////////////////////////////////////////////////////////


// %%COMMENT
// TODO: Adjust this to new order of basis vectors rep!!!

#include <string.h>
#include "mm_op15.h"   
   




static const uint_mmv_t TABLE_PERM64_LOW[] = {
   // %%TABLE TABLE_PERM64_XY_LOW, uint{INT_BITS}
0x0000000000000000ULL,0x0000000000000000ULL,
0x0000000000000000ULL,0x0000000000000000ULL,
0xff0000ff00ffff00ULL,0x00ffff00ff0000ffULL,
0x00ffff00ff0000ffULL,0xff0000ff00ffff00ULL,
0xf0f00f0f0f0ff0f0ULL,0x0f0ff0f0f0f00f0fULL,
0x0f0ff0f0f0f00f0fULL,0xf0f00f0f0f0ff0f0ULL,
0x0ff00ff00ff00ff0ULL,0x0ff00ff00ff00ff0ULL,
0x0ff00ff00ff00ff0ULL,0x0ff00ff00ff00ff0ULL,
0xf00ff00f0ff00ff0ULL,0x0ff00ff0f00ff00fULL,
0x0ff00ff0f00ff00fULL,0xf00ff00f0ff00ff0ULL,
0x0f0ff0f00f0ff0f0ULL,0x0f0ff0f00f0ff0f0ULL,
0x0f0ff0f00f0ff0f0ULL,0x0f0ff0f00f0ff0f0ULL,
0x00ffff0000ffff00ULL,0x00ffff0000ffff00ULL,
0x00ffff0000ffff00ULL,0x00ffff0000ffff00ULL,
0xffffffff00000000ULL,0x00000000ffffffffULL,
0x00000000ffffffffULL,0xffffffff00000000ULL,
0xf00f0ff0f00f0ff0ULL,0x0ff0f00f0ff0f00fULL,
0x0ff0f00f0ff0f00fULL,0xf00f0ff0f00f0ff0ULL,
0x0f0f0f0ff0f0f0f0ULL,0x0f0f0f0ff0f0f0f0ULL,
0x0f0f0f0ff0f0f0f0ULL,0x0f0f0f0ff0f0f0f0ULL,
0x00ff00ffff00ff00ULL,0x00ff00ffff00ff00ULL,
0x00ff00ffff00ff00ULL,0x00ff00ffff00ff00ULL,
0xffff0000ffff0000ULL,0x0000ffff0000ffffULL,
0x0000ffff0000ffffULL,0xffff0000ffff0000ULL,
0x0000ffffffff0000ULL,0x0000ffffffff0000ULL,
0x0000ffffffff0000ULL,0x0000ffffffff0000ULL,
0xff00ff00ff00ff00ULL,0x00ff00ff00ff00ffULL,
0x00ff00ff00ff00ffULL,0xff00ff00ff00ff00ULL,
0xf0f0f0f0f0f0f0f0ULL,0x0f0f0f0f0f0f0f0fULL,
0x0f0f0f0f0f0f0f0fULL,0xf0f0f0f0f0f0f0f0ULL,
0x0ff0f00ff00f0ff0ULL,0x0ff0f00ff00f0ff0ULL,
0x0ff0f00ff00f0ff0ULL,0x0ff0f00ff00f0ff0ULL
}; 


static const uint_mmv_t TABLE_PERM64_HIGH[] = {
   // %%TABLE TABLE_PERM64_XY_HIGH, uint{INT_BITS}
0x0000000000000000ULL,0x0000000000000000ULL,
0x0000000000000000ULL,0x0000000000000000ULL,
0x0ff0f00ff00f0ff0ULL,0x0ff0f00ff00f0ff0ULL,
0xf00f0ff00ff0f00fULL,0xf00f0ff00ff0f00fULL,
0x0ff0f00ff00f0ff0ULL,0xf00f0ff00ff0f00fULL,
0x0ff0f00ff00f0ff0ULL,0xf00f0ff00ff0f00fULL,
0x0000000000000000ULL,0xffffffffffffffffULL,
0xffffffffffffffffULL,0x0000000000000000ULL,
0xffffffffffffffffULL,0xffffffffffffffffULL,
0xffffffffffffffffULL,0xffffffffffffffffULL,
0xf00f0ff00ff0f00fULL,0xf00f0ff00ff0f00fULL,
0x0ff0f00ff00f0ff0ULL,0x0ff0f00ff00f0ff0ULL,
0xf00f0ff00ff0f00fULL,0x0ff0f00ff00f0ff0ULL,
0xf00f0ff00ff0f00fULL,0x0ff0f00ff00f0ff0ULL,
0xffffffffffffffffULL,0x0000000000000000ULL,
0x0000000000000000ULL,0xffffffffffffffffULL,
0x000f0fff0ffffff0ULL,0xf000000f000f0fffULL,
0xf000000f000f0fffULL,0xfff0f000f000000fULL,
0x0ffffff0fff0f000ULL,0xfff0f000f000000fULL,
0x000f0fff0ffffff0ULL,0x0ffffff0fff0f000ULL,
0x0ffffff0fff0f000ULL,0x000f0fff0ffffff0ULL,
0xfff0f000f000000fULL,0x0ffffff0fff0f000ULL,
0x000f0fff0ffffff0ULL,0x0ffffff0fff0f000ULL,
0x0ffffff0fff0f000ULL,0xfff0f000f000000fULL,
0xfff0f000f000000fULL,0x0ffffff0fff0f000ULL,
0x0ffffff0fff0f000ULL,0x000f0fff0ffffff0ULL,
0xf000000f000f0fffULL,0x000f0fff0ffffff0ULL,
0xfff0f000f000000fULL,0xf000000f000f0fffULL,
0xf000000f000f0fffULL,0xfff0f000f000000fULL,
0x000f0fff0ffffff0ULL,0xf000000f000f0fffULL,
0xfff0f000f000000fULL,0xf000000f000f0fffULL,
0xf000000f000f0fffULL,0x000f0fff0ffffff0ULL
}; 


static const uint32_t TABLE24_START[4] = {
   MM_OP15_OFS_X, MM_OP15_OFS_Z, MM_OP15_OFS_Y, MM_OP15_OFS_A
};




// %%EXPORT
void mm_op15_do_xy(uint_mmv_t *v_in, mm_sub_op_xy_type *p_op, uint_mmv_t *v_out)
{
    uint_fast32_t i;

    for (i = 0; i < MM_OP15_OFS_E; ++i) v_out[i] = 0;
    
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
            uint_mmv_t a_sign[2][2];
            uint_mmv_t d_xor = p_op->lin_d[i];
            uint8_t *p_sign = p_op->sign_XYZ;
    
            for (i1 = 0; i1 < 2; ++i1) {
                uint_mmv_t x = p_op->lin_i[i] >> (i1 << 4); 
                // %%MMV_UINT_SPREAD x, x
                // Spread bits 0,...,15 of x to the (4-bit long) fields
                // of x. A field of x is set to 0xf if its 
                // corresponding bit in input x is one and to 0 otherwise.
                x = (x & 0xffULL) + ((x & 0xff00ULL) << 24);
                x = (x & 0xf0000000fULL) 
                    +  ((x & 0xf0000000f0ULL) << 12);
                x = (x & 0x3000300030003ULL) 
                    +  ((x & 0xc000c000c000cULL) << 6);
                x = (x & 0x101010101010101ULL) 
                    +  ((x & 0x202020202020202ULL) << 3);
                x = (((x) << 4) - (x));
                // Bit spreading done.
                a_sign[0][i1] = x;
                a_sign[1][i1] = x ^ 0xffffffffffffffffULL;
            }
            a_sign[1][1] &= 0xffffffffULL;
         
            for (i1 = 0; i1 < 2048; ++i1) {
                uint_mmv_t *ps = p_src + ((i1 ^ d_xor) << 1);
                uint_fast8_t sign = (p_sign[i1] >> i) & 1;
                // %%FOR j in range(V24_INTS_USED)
                p_dest[0] = ps[0] ^ a_sign[sign][0];
                p_dest[1] = ps[1] ^ a_sign[sign][1];
                // %%END FOR
                p_dest +=  2;      
            }
        }    
    }    

    // Step 2: do rows with 64 entries, tag T // TODO: comment properly!!!!
    {
        uint_mmv_t *p_src = v_in + MM_OP15_OFS_T;
        uint_mmv_t *p_dest = v_out + MM_OP15_OFS_T;
        uint16_t* p_T =  p_op->s_T;
        for (i = 0; i < 759; ++i) {
            uint_fast16_t ofs_l = *p_T;
            uint_fast16_t ofs_h = (ofs_l & 63) >> 4;
            const uint_mmv_t *ps_h = TABLE_PERM64_HIGH +
                ((ofs_l & 0xf000) >> 10);
            const uint_mmv_t *ps_l = TABLE_PERM64_LOW + 
                ((ofs_l & 0xf00) >> 6);
            ofs_l = (ofs_l << 2) & 0x3fULL;
            // %%FOR j in range(V64_INTS)
            p_dest[0] =  ps_h[0] ^ ps_l[0] ^
               (((p_src[0 ^ ofs_h] >> (0 ^ ofs_l)) & 15) << 0) ^
               (((p_src[0 ^ ofs_h] >> (4 ^ ofs_l)) & 15) << 4) ^
               (((p_src[0 ^ ofs_h] >> (8 ^ ofs_l)) & 15) << 8) ^
               (((p_src[0 ^ ofs_h] >> (12 ^ ofs_l)) & 15) << 12) ^
               (((p_src[0 ^ ofs_h] >> (16 ^ ofs_l)) & 15) << 16) ^
               (((p_src[0 ^ ofs_h] >> (20 ^ ofs_l)) & 15) << 20) ^
               (((p_src[0 ^ ofs_h] >> (24 ^ ofs_l)) & 15) << 24) ^
               (((p_src[0 ^ ofs_h] >> (28 ^ ofs_l)) & 15) << 28) ^
               (((p_src[0 ^ ofs_h] >> (32 ^ ofs_l)) & 15) << 32) ^
               (((p_src[0 ^ ofs_h] >> (36 ^ ofs_l)) & 15) << 36) ^
               (((p_src[0 ^ ofs_h] >> (40 ^ ofs_l)) & 15) << 40) ^
               (((p_src[0 ^ ofs_h] >> (44 ^ ofs_l)) & 15) << 44) ^
               (((p_src[0 ^ ofs_h] >> (48 ^ ofs_l)) & 15) << 48) ^
               (((p_src[0 ^ ofs_h] >> (52 ^ ofs_l)) & 15) << 52) ^
               (((p_src[0 ^ ofs_h] >> (56 ^ ofs_l)) & 15) << 56) ^
               (((p_src[0 ^ ofs_h] >> (60 ^ ofs_l)) & 15) << 60);
            p_dest[1] =  ps_h[1] ^ ps_l[1] ^
               (((p_src[1 ^ ofs_h] >> (0 ^ ofs_l)) & 15) << 0) ^
               (((p_src[1 ^ ofs_h] >> (4 ^ ofs_l)) & 15) << 4) ^
               (((p_src[1 ^ ofs_h] >> (8 ^ ofs_l)) & 15) << 8) ^
               (((p_src[1 ^ ofs_h] >> (12 ^ ofs_l)) & 15) << 12) ^
               (((p_src[1 ^ ofs_h] >> (16 ^ ofs_l)) & 15) << 16) ^
               (((p_src[1 ^ ofs_h] >> (20 ^ ofs_l)) & 15) << 20) ^
               (((p_src[1 ^ ofs_h] >> (24 ^ ofs_l)) & 15) << 24) ^
               (((p_src[1 ^ ofs_h] >> (28 ^ ofs_l)) & 15) << 28) ^
               (((p_src[1 ^ ofs_h] >> (32 ^ ofs_l)) & 15) << 32) ^
               (((p_src[1 ^ ofs_h] >> (36 ^ ofs_l)) & 15) << 36) ^
               (((p_src[1 ^ ofs_h] >> (40 ^ ofs_l)) & 15) << 40) ^
               (((p_src[1 ^ ofs_h] >> (44 ^ ofs_l)) & 15) << 44) ^
               (((p_src[1 ^ ofs_h] >> (48 ^ ofs_l)) & 15) << 48) ^
               (((p_src[1 ^ ofs_h] >> (52 ^ ofs_l)) & 15) << 52) ^
               (((p_src[1 ^ ofs_h] >> (56 ^ ofs_l)) & 15) << 56) ^
               (((p_src[1 ^ ofs_h] >> (60 ^ ofs_l)) & 15) << 60);
            p_dest[2] =  ps_h[2] ^ ps_l[2] ^
               (((p_src[2 ^ ofs_h] >> (0 ^ ofs_l)) & 15) << 0) ^
               (((p_src[2 ^ ofs_h] >> (4 ^ ofs_l)) & 15) << 4) ^
               (((p_src[2 ^ ofs_h] >> (8 ^ ofs_l)) & 15) << 8) ^
               (((p_src[2 ^ ofs_h] >> (12 ^ ofs_l)) & 15) << 12) ^
               (((p_src[2 ^ ofs_h] >> (16 ^ ofs_l)) & 15) << 16) ^
               (((p_src[2 ^ ofs_h] >> (20 ^ ofs_l)) & 15) << 20) ^
               (((p_src[2 ^ ofs_h] >> (24 ^ ofs_l)) & 15) << 24) ^
               (((p_src[2 ^ ofs_h] >> (28 ^ ofs_l)) & 15) << 28) ^
               (((p_src[2 ^ ofs_h] >> (32 ^ ofs_l)) & 15) << 32) ^
               (((p_src[2 ^ ofs_h] >> (36 ^ ofs_l)) & 15) << 36) ^
               (((p_src[2 ^ ofs_h] >> (40 ^ ofs_l)) & 15) << 40) ^
               (((p_src[2 ^ ofs_h] >> (44 ^ ofs_l)) & 15) << 44) ^
               (((p_src[2 ^ ofs_h] >> (48 ^ ofs_l)) & 15) << 48) ^
               (((p_src[2 ^ ofs_h] >> (52 ^ ofs_l)) & 15) << 52) ^
               (((p_src[2 ^ ofs_h] >> (56 ^ ofs_l)) & 15) << 56) ^
               (((p_src[2 ^ ofs_h] >> (60 ^ ofs_l)) & 15) << 60);
            p_dest[3] =  ps_h[3] ^ ps_l[3] ^
               (((p_src[3 ^ ofs_h] >> (0 ^ ofs_l)) & 15) << 0) ^
               (((p_src[3 ^ ofs_h] >> (4 ^ ofs_l)) & 15) << 4) ^
               (((p_src[3 ^ ofs_h] >> (8 ^ ofs_l)) & 15) << 8) ^
               (((p_src[3 ^ ofs_h] >> (12 ^ ofs_l)) & 15) << 12) ^
               (((p_src[3 ^ ofs_h] >> (16 ^ ofs_l)) & 15) << 16) ^
               (((p_src[3 ^ ofs_h] >> (20 ^ ofs_l)) & 15) << 20) ^
               (((p_src[3 ^ ofs_h] >> (24 ^ ofs_l)) & 15) << 24) ^
               (((p_src[3 ^ ofs_h] >> (28 ^ ofs_l)) & 15) << 28) ^
               (((p_src[3 ^ ofs_h] >> (32 ^ ofs_l)) & 15) << 32) ^
               (((p_src[3 ^ ofs_h] >> (36 ^ ofs_l)) & 15) << 36) ^
               (((p_src[3 ^ ofs_h] >> (40 ^ ofs_l)) & 15) << 40) ^
               (((p_src[3 ^ ofs_h] >> (44 ^ ofs_l)) & 15) << 44) ^
               (((p_src[3 ^ ofs_h] >> (48 ^ ofs_l)) & 15) << 48) ^
               (((p_src[3 ^ ofs_h] >> (52 ^ ofs_l)) & 15) << 52) ^
               (((p_src[3 ^ ofs_h] >> (56 ^ ofs_l)) & 15) << 56) ^
               (((p_src[3 ^ ofs_h] >> (60 ^ ofs_l)) & 15) << 60);
            // %%END FOR
            p_src += 4; 
            p_dest += 4; 
            ++p_T;
        }
    }


    // Step 3: do rows with 24 entries, tags A, B, C // TODO: comment properly!!!!
    {
        uint_mmv_t mask[8];
        uint_mmv_t neg_mask[2];
        uint_mmv_t f = p_op->f_i, ef = p_op->ef_i, eps;
        for (i = 0; i < 2; ++i) {
             mask[i] = f >> (i << 4);
             mask[i + 2] = ef >> (i << 4);
        }
        neg_mask[0] = 0xffffffffffffffffULL;
        neg_mask[1] = 0xffffffffULL;
        for (i = 0; i < 4; ++i) {
            uint_mmv_t x = mask[i];
            // %%MMV_UINT_SPREAD x, x
            // Spread bits 0,...,15 of x to the (4-bit long) fields
            // of x. A field of x is set to 0xf if its 
            // corresponding bit in input x is one and to 0 otherwise.
            x = (x & 0xffULL) + ((x & 0xff00ULL) << 24);
            x = (x & 0xf0000000fULL) 
                +  ((x & 0xf0000000f0ULL) << 12);
            x = (x & 0x3000300030003ULL) 
                +  ((x & 0xc000c000c000cULL) << 6);
            x = (x & 0x101010101010101ULL) 
                +  ((x & 0x202020202020202ULL) << 3);
            x = (((x) << 4) - (x));
            // Bit spreading done.
            mask[i] = x = x & neg_mask[i & 1];
            mask[i + 4] = x ^ neg_mask[i & 1];
        }

        f <<= 2;
        ef <<= 2;
        eps = -(p_op->eps & 0x1ULL);
        for (i = 0; i < 48; i += 2) {
            uint_mmv_t t, t1, t2;
            // %%FOR j in range(V24_INTS_USED)
            // process uint_mmv_t 0 of row i/2 for tags A, B, C
            t1 = v_in[i + MM_OP15_OFS_A + 0]; 
            t = mask[0 + (f & 0x4ULL)];
            v_out[i + MM_OP15_OFS_A + 0] = t1 ^ t; 
            t1 = v_in[i + MM_OP15_OFS_B + 0]; 
            t2 = v_in[i + MM_OP15_OFS_C + 0];
            t &= (t1 ^ t2);
            t ^= mask[2 + (ef & 0x4ULL)];
            v_out[i + MM_OP15_OFS_B + 0] = t1 ^ t;
            t2 ^= eps & neg_mask[0];
            v_out[i + MM_OP15_OFS_C + 0] = t2 ^ t;
            // process uint_mmv_t 1 of row i/2 for tags A, B, C
            t1 = v_in[i + MM_OP15_OFS_A + 1]; 
            t = mask[1 + (f & 0x4ULL)];
            v_out[i + MM_OP15_OFS_A + 1] = t1 ^ t; 
            t1 = v_in[i + MM_OP15_OFS_B + 1]; 
            t2 = v_in[i + MM_OP15_OFS_C + 1];
            t &= (t1 ^ t2);
            t ^= mask[3 + (ef & 0x4ULL)];
            v_out[i + MM_OP15_OFS_B + 1] = t1 ^ t;
            t2 ^= eps & neg_mask[1];
            v_out[i + MM_OP15_OFS_C + 1] = t2 ^ t;
            // %%END FOR
            f >>= 1; ef >>= 1;      
        }
        //yet to be done!!!!
    }



    // If eps is odd: 
    //    negate entries X_d,i with scalar product <d,i> = 1
    if (p_op->eps & 1) mm15_neg_scalprod_d_i(v_out + MM_OP15_OFS_X); 
} 





// %%EXPORT p
void mm_op15_xy(uint_mmv_t *v_in, uint32_t f, uint32_t e, uint32_t eps, uint_mmv_t *v_out)
{
    mm_sub_op_xy_type s_op;
    mm_sub_prep_xy(f, e, eps, &s_op);
    mm_op15_do_xy(v_in, &s_op, v_out);
}


