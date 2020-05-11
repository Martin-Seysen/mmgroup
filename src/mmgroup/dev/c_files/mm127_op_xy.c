/////////////////////////////////////////////////////////////////////////////
// This C file has been created automatically. Do not edit!!!
/////////////////////////////////////////////////////////////////////////////


// %%COMMENT
// TODO: Adjust this to new order of basis vectors rep!!!

#include <string.h>
#include "mm_op127.h"   
   




static const uint_mmv_t TABLE_PERM64_LOW[] = {
   // %%TABLE TABLE_PERM64_XY_LOW, uint{INT_BITS}
0x0000000000000000ULL,0x0000000000000000ULL,
0x0000000000000000ULL,0x0000000000000000ULL,
0x0000000000000000ULL,0x0000000000000000ULL,
0x0000000000000000ULL,0x0000000000000000ULL,
0x00007f7f7f7f0000ULL,0x7f7f000000007f7fULL,
0x7f7f000000007f7fULL,0x00007f7f7f7f0000ULL,
0x7f7f000000007f7fULL,0x00007f7f7f7f0000ULL,
0x00007f7f7f7f0000ULL,0x7f7f000000007f7fULL,
0x007f007f7f007f00ULL,0x7f007f00007f007fULL,
0x7f007f00007f007fULL,0x007f007f7f007f00ULL,
0x7f007f00007f007fULL,0x007f007f7f007f00ULL,
0x007f007f7f007f00ULL,0x7f007f00007f007fULL,
0x007f7f00007f7f00ULL,0x007f7f00007f7f00ULL,
0x007f7f00007f7f00ULL,0x007f7f00007f7f00ULL,
0x007f7f00007f7f00ULL,0x007f7f00007f7f00ULL,
0x007f7f00007f7f00ULL,0x007f7f00007f7f00ULL,
0x007f7f00007f7f00ULL,0x7f00007f7f00007fULL,
0x7f00007f7f00007fULL,0x007f7f00007f7f00ULL,
0x7f00007f7f00007fULL,0x007f7f00007f7f00ULL,
0x007f7f00007f7f00ULL,0x7f00007f7f00007fULL,
0x007f007f7f007f00ULL,0x007f007f7f007f00ULL,
0x007f007f7f007f00ULL,0x007f007f7f007f00ULL,
0x007f007f7f007f00ULL,0x007f007f7f007f00ULL,
0x007f007f7f007f00ULL,0x007f007f7f007f00ULL,
0x00007f7f7f7f0000ULL,0x00007f7f7f7f0000ULL,
0x00007f7f7f7f0000ULL,0x00007f7f7f7f0000ULL,
0x00007f7f7f7f0000ULL,0x00007f7f7f7f0000ULL,
0x00007f7f7f7f0000ULL,0x00007f7f7f7f0000ULL,
0x0000000000000000ULL,0x7f7f7f7f7f7f7f7fULL,
0x7f7f7f7f7f7f7f7fULL,0x0000000000000000ULL,
0x7f7f7f7f7f7f7f7fULL,0x0000000000000000ULL,
0x0000000000000000ULL,0x7f7f7f7f7f7f7f7fULL,
0x7f00007f007f7f00ULL,0x7f00007f007f7f00ULL,
0x007f7f007f00007fULL,0x007f7f007f00007fULL,
0x007f7f007f00007fULL,0x007f7f007f00007fULL,
0x7f00007f007f7f00ULL,0x7f00007f007f7f00ULL,
0x7f007f007f007f00ULL,0x007f007f007f007fULL,
0x7f007f007f007f00ULL,0x007f007f007f007fULL,
0x7f007f007f007f00ULL,0x007f007f007f007fULL,
0x7f007f007f007f00ULL,0x007f007f007f007fULL,
0x7f7f00007f7f0000ULL,0x00007f7f00007f7fULL,
0x7f7f00007f7f0000ULL,0x00007f7f00007f7fULL,
0x7f7f00007f7f0000ULL,0x00007f7f00007f7fULL,
0x7f7f00007f7f0000ULL,0x00007f7f00007f7fULL,
0x7f7f7f7f00000000ULL,0x7f7f7f7f00000000ULL,
0x000000007f7f7f7fULL,0x000000007f7f7f7fULL,
0x000000007f7f7f7fULL,0x000000007f7f7f7fULL,
0x7f7f7f7f00000000ULL,0x7f7f7f7f00000000ULL,
0x7f7f7f7f00000000ULL,0x000000007f7f7f7fULL,
0x7f7f7f7f00000000ULL,0x000000007f7f7f7fULL,
0x7f7f7f7f00000000ULL,0x000000007f7f7f7fULL,
0x7f7f7f7f00000000ULL,0x000000007f7f7f7fULL,
0x7f7f00007f7f0000ULL,0x7f7f00007f7f0000ULL,
0x00007f7f00007f7fULL,0x00007f7f00007f7fULL,
0x00007f7f00007f7fULL,0x00007f7f00007f7fULL,
0x7f7f00007f7f0000ULL,0x7f7f00007f7f0000ULL,
0x7f007f007f007f00ULL,0x7f007f007f007f00ULL,
0x007f007f007f007fULL,0x007f007f007f007fULL,
0x007f007f007f007fULL,0x007f007f007f007fULL,
0x7f007f007f007f00ULL,0x7f007f007f007f00ULL,
0x7f00007f007f7f00ULL,0x007f7f007f00007fULL,
0x7f00007f007f7f00ULL,0x007f7f007f00007fULL,
0x7f00007f007f7f00ULL,0x007f7f007f00007fULL,
0x7f00007f007f7f00ULL,0x007f7f007f00007fULL
}; 


static const uint_mmv_t TABLE_PERM64_HIGH[] = {
   // %%TABLE TABLE_PERM64_XY_HIGH, uint{INT_BITS}
0x0000000000000000ULL,0x0000000000000000ULL,
0x0000000000000000ULL,0x0000000000000000ULL,
0x0000000000000000ULL,0x0000000000000000ULL,
0x0000000000000000ULL,0x0000000000000000ULL,
0x7f00007f007f7f00ULL,0x007f7f007f00007fULL,
0x7f00007f007f7f00ULL,0x007f7f007f00007fULL,
0x007f7f007f00007fULL,0x7f00007f007f7f00ULL,
0x007f7f007f00007fULL,0x7f00007f007f7f00ULL,
0x7f00007f007f7f00ULL,0x007f7f007f00007fULL,
0x007f7f007f00007fULL,0x7f00007f007f7f00ULL,
0x7f00007f007f7f00ULL,0x007f7f007f00007fULL,
0x007f7f007f00007fULL,0x7f00007f007f7f00ULL,
0x0000000000000000ULL,0x0000000000000000ULL,
0x7f7f7f7f7f7f7f7fULL,0x7f7f7f7f7f7f7f7fULL,
0x7f7f7f7f7f7f7f7fULL,0x7f7f7f7f7f7f7f7fULL,
0x0000000000000000ULL,0x0000000000000000ULL,
0x7f7f7f7f7f7f7f7fULL,0x7f7f7f7f7f7f7f7fULL,
0x7f7f7f7f7f7f7f7fULL,0x7f7f7f7f7f7f7f7fULL,
0x7f7f7f7f7f7f7f7fULL,0x7f7f7f7f7f7f7f7fULL,
0x7f7f7f7f7f7f7f7fULL,0x7f7f7f7f7f7f7f7fULL,
0x007f7f007f00007fULL,0x7f00007f007f7f00ULL,
0x007f7f007f00007fULL,0x7f00007f007f7f00ULL,
0x7f00007f007f7f00ULL,0x007f7f007f00007fULL,
0x7f00007f007f7f00ULL,0x007f7f007f00007fULL,
0x007f7f007f00007fULL,0x7f00007f007f7f00ULL,
0x7f00007f007f7f00ULL,0x007f7f007f00007fULL,
0x007f7f007f00007fULL,0x7f00007f007f7f00ULL,
0x7f00007f007f7f00ULL,0x007f7f007f00007fULL,
0x7f7f7f7f7f7f7f7fULL,0x7f7f7f7f7f7f7f7fULL,
0x0000000000000000ULL,0x0000000000000000ULL,
0x0000000000000000ULL,0x0000000000000000ULL,
0x7f7f7f7f7f7f7f7fULL,0x7f7f7f7f7f7f7f7fULL,
0x007f7f7f7f7f7f00ULL,0x0000007f007f7f7fULL,
0x0000007f007f7f7fULL,0x7f0000000000007fULL,
0x0000007f007f7f7fULL,0x7f0000000000007fULL,
0x7f0000000000007fULL,0x7f7f7f007f000000ULL,
0x7f7f7f007f000000ULL,0x007f7f7f7f7f7f00ULL,
0x7f0000000000007fULL,0x7f7f7f007f000000ULL,
0x007f7f7f7f7f7f00ULL,0x0000007f007f7f7fULL,
0x7f7f7f007f000000ULL,0x007f7f7f7f7f7f00ULL,
0x7f7f7f007f000000ULL,0x007f7f7f7f7f7f00ULL,
0x007f7f7f7f7f7f00ULL,0x0000007f007f7f7fULL,
0x7f0000000000007fULL,0x7f7f7f007f000000ULL,
0x7f7f7f007f000000ULL,0x007f7f7f7f7f7f00ULL,
0x007f7f7f7f7f7f00ULL,0x0000007f007f7f7fULL,
0x7f7f7f007f000000ULL,0x007f7f7f7f7f7f00ULL,
0x7f7f7f007f000000ULL,0x007f7f7f7f7f7f00ULL,
0x7f0000000000007fULL,0x7f7f7f007f000000ULL,
0x7f0000000000007fULL,0x7f7f7f007f000000ULL,
0x7f7f7f007f000000ULL,0x007f7f7f7f7f7f00ULL,
0x7f7f7f007f000000ULL,0x007f7f7f7f7f7f00ULL,
0x007f7f7f7f7f7f00ULL,0x0000007f007f7f7fULL,
0x0000007f007f7f7fULL,0x7f0000000000007fULL,
0x007f7f7f7f7f7f00ULL,0x0000007f007f7f7fULL,
0x7f0000000000007fULL,0x7f7f7f007f000000ULL,
0x0000007f007f7f7fULL,0x7f0000000000007fULL,
0x0000007f007f7f7fULL,0x7f0000000000007fULL,
0x7f0000000000007fULL,0x7f7f7f007f000000ULL,
0x007f7f7f7f7f7f00ULL,0x0000007f007f7f7fULL,
0x0000007f007f7f7fULL,0x7f0000000000007fULL,
0x7f0000000000007fULL,0x7f7f7f007f000000ULL,
0x0000007f007f7f7fULL,0x7f0000000000007fULL,
0x0000007f007f7f7fULL,0x7f0000000000007fULL,
0x007f7f7f7f7f7f00ULL,0x0000007f007f7f7fULL
}; 


static const uint32_t TABLE24_START[4] = {
   MM_OP127_OFS_X, MM_OP127_OFS_Z, MM_OP127_OFS_Y, MM_OP127_OFS_A
};




// %%EXPORT
void mm_op127_do_xy(uint_mmv_t *v_in, mm_sub_op_xy_type *p_op, uint_mmv_t *v_out)
{
    uint_fast32_t i;

    for (i = 0; i < MM_OP127_OFS_E; ++i) v_out[i] = 0;
    
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
            uint_mmv_t a_sign[2][4];
            uint_mmv_t d_xor = p_op->lin_d[i];
            uint8_t *p_sign = p_op->sign_XYZ;
    
            for (i1 = 0; i1 < 3; ++i1) {
                uint_mmv_t x = p_op->lin_i[i] >> (i1 << 3); 
                // %%MMV_UINT_SPREAD x, x
                // Spread bits 0,...,7 of x to the (8-bit long) fields
                // of x. A field of x is set to 0x7f if its 
                // corresponding bit in input x is one and to 0 otherwise.
                x = (x & 0xfULL) + ((x & 0xf0ULL) << 28);
                x = (x & 0x300000003ULL) 
                    +  ((x & 0xc0000000cULL) << 14);
                x = (x & 0x1000100010001ULL) 
                    +  ((x & 0x2000200020002ULL) << 7);
                x = (((x) << 7) - (x));
                // Bit spreading done.
                a_sign[0][i1] = x;
                a_sign[1][i1] = x ^ 0x7f7f7f7f7f7f7f7fULL;
            }
         
            for (i1 = 0; i1 < 2048; ++i1) {
                uint_mmv_t *ps = p_src + ((i1 ^ d_xor) << 2);
                uint_fast8_t sign = (p_sign[i1] >> i) & 1;
                // %%FOR j in range(V24_INTS_USED)
                p_dest[0] = ps[0] ^ a_sign[sign][0];
                p_dest[1] = ps[1] ^ a_sign[sign][1];
                p_dest[2] = ps[2] ^ a_sign[sign][2];
                // %%END FOR
                p_dest[3] = 0;
                p_dest +=  4;      
            }
        }    
    }    

    // Step 2: do rows with 64 entries, tag T // TODO: comment properly!!!!
    {
        uint_mmv_t *p_src = v_in + MM_OP127_OFS_T;
        uint_mmv_t *p_dest = v_out + MM_OP127_OFS_T;
        uint16_t* p_T =  p_op->s_T;
        for (i = 0; i < 759; ++i) {
            uint_fast16_t ofs_l = *p_T;
            uint_fast16_t ofs_h = (ofs_l & 63) >> 3;
            const uint_mmv_t *ps_h = TABLE_PERM64_HIGH +
                ((ofs_l & 0xf000) >> 9);
            const uint_mmv_t *ps_l = TABLE_PERM64_LOW + 
                ((ofs_l & 0xf00) >> 5);
            ofs_l = (ofs_l << 3) & 0x3fULL;
            // %%FOR j in range(V64_INTS)
            p_dest[0] =  ps_h[0] ^ ps_l[0] ^
               (((p_src[0 ^ ofs_h] >> (0 ^ ofs_l)) & 127) << 0) ^
               (((p_src[0 ^ ofs_h] >> (8 ^ ofs_l)) & 127) << 8) ^
               (((p_src[0 ^ ofs_h] >> (16 ^ ofs_l)) & 127) << 16) ^
               (((p_src[0 ^ ofs_h] >> (24 ^ ofs_l)) & 127) << 24) ^
               (((p_src[0 ^ ofs_h] >> (32 ^ ofs_l)) & 127) << 32) ^
               (((p_src[0 ^ ofs_h] >> (40 ^ ofs_l)) & 127) << 40) ^
               (((p_src[0 ^ ofs_h] >> (48 ^ ofs_l)) & 127) << 48) ^
               (((p_src[0 ^ ofs_h] >> (56 ^ ofs_l)) & 127) << 56);
            p_dest[1] =  ps_h[1] ^ ps_l[1] ^
               (((p_src[1 ^ ofs_h] >> (0 ^ ofs_l)) & 127) << 0) ^
               (((p_src[1 ^ ofs_h] >> (8 ^ ofs_l)) & 127) << 8) ^
               (((p_src[1 ^ ofs_h] >> (16 ^ ofs_l)) & 127) << 16) ^
               (((p_src[1 ^ ofs_h] >> (24 ^ ofs_l)) & 127) << 24) ^
               (((p_src[1 ^ ofs_h] >> (32 ^ ofs_l)) & 127) << 32) ^
               (((p_src[1 ^ ofs_h] >> (40 ^ ofs_l)) & 127) << 40) ^
               (((p_src[1 ^ ofs_h] >> (48 ^ ofs_l)) & 127) << 48) ^
               (((p_src[1 ^ ofs_h] >> (56 ^ ofs_l)) & 127) << 56);
            p_dest[2] =  ps_h[2] ^ ps_l[2] ^
               (((p_src[2 ^ ofs_h] >> (0 ^ ofs_l)) & 127) << 0) ^
               (((p_src[2 ^ ofs_h] >> (8 ^ ofs_l)) & 127) << 8) ^
               (((p_src[2 ^ ofs_h] >> (16 ^ ofs_l)) & 127) << 16) ^
               (((p_src[2 ^ ofs_h] >> (24 ^ ofs_l)) & 127) << 24) ^
               (((p_src[2 ^ ofs_h] >> (32 ^ ofs_l)) & 127) << 32) ^
               (((p_src[2 ^ ofs_h] >> (40 ^ ofs_l)) & 127) << 40) ^
               (((p_src[2 ^ ofs_h] >> (48 ^ ofs_l)) & 127) << 48) ^
               (((p_src[2 ^ ofs_h] >> (56 ^ ofs_l)) & 127) << 56);
            p_dest[3] =  ps_h[3] ^ ps_l[3] ^
               (((p_src[3 ^ ofs_h] >> (0 ^ ofs_l)) & 127) << 0) ^
               (((p_src[3 ^ ofs_h] >> (8 ^ ofs_l)) & 127) << 8) ^
               (((p_src[3 ^ ofs_h] >> (16 ^ ofs_l)) & 127) << 16) ^
               (((p_src[3 ^ ofs_h] >> (24 ^ ofs_l)) & 127) << 24) ^
               (((p_src[3 ^ ofs_h] >> (32 ^ ofs_l)) & 127) << 32) ^
               (((p_src[3 ^ ofs_h] >> (40 ^ ofs_l)) & 127) << 40) ^
               (((p_src[3 ^ ofs_h] >> (48 ^ ofs_l)) & 127) << 48) ^
               (((p_src[3 ^ ofs_h] >> (56 ^ ofs_l)) & 127) << 56);
            p_dest[4] =  ps_h[4] ^ ps_l[4] ^
               (((p_src[4 ^ ofs_h] >> (0 ^ ofs_l)) & 127) << 0) ^
               (((p_src[4 ^ ofs_h] >> (8 ^ ofs_l)) & 127) << 8) ^
               (((p_src[4 ^ ofs_h] >> (16 ^ ofs_l)) & 127) << 16) ^
               (((p_src[4 ^ ofs_h] >> (24 ^ ofs_l)) & 127) << 24) ^
               (((p_src[4 ^ ofs_h] >> (32 ^ ofs_l)) & 127) << 32) ^
               (((p_src[4 ^ ofs_h] >> (40 ^ ofs_l)) & 127) << 40) ^
               (((p_src[4 ^ ofs_h] >> (48 ^ ofs_l)) & 127) << 48) ^
               (((p_src[4 ^ ofs_h] >> (56 ^ ofs_l)) & 127) << 56);
            p_dest[5] =  ps_h[5] ^ ps_l[5] ^
               (((p_src[5 ^ ofs_h] >> (0 ^ ofs_l)) & 127) << 0) ^
               (((p_src[5 ^ ofs_h] >> (8 ^ ofs_l)) & 127) << 8) ^
               (((p_src[5 ^ ofs_h] >> (16 ^ ofs_l)) & 127) << 16) ^
               (((p_src[5 ^ ofs_h] >> (24 ^ ofs_l)) & 127) << 24) ^
               (((p_src[5 ^ ofs_h] >> (32 ^ ofs_l)) & 127) << 32) ^
               (((p_src[5 ^ ofs_h] >> (40 ^ ofs_l)) & 127) << 40) ^
               (((p_src[5 ^ ofs_h] >> (48 ^ ofs_l)) & 127) << 48) ^
               (((p_src[5 ^ ofs_h] >> (56 ^ ofs_l)) & 127) << 56);
            p_dest[6] =  ps_h[6] ^ ps_l[6] ^
               (((p_src[6 ^ ofs_h] >> (0 ^ ofs_l)) & 127) << 0) ^
               (((p_src[6 ^ ofs_h] >> (8 ^ ofs_l)) & 127) << 8) ^
               (((p_src[6 ^ ofs_h] >> (16 ^ ofs_l)) & 127) << 16) ^
               (((p_src[6 ^ ofs_h] >> (24 ^ ofs_l)) & 127) << 24) ^
               (((p_src[6 ^ ofs_h] >> (32 ^ ofs_l)) & 127) << 32) ^
               (((p_src[6 ^ ofs_h] >> (40 ^ ofs_l)) & 127) << 40) ^
               (((p_src[6 ^ ofs_h] >> (48 ^ ofs_l)) & 127) << 48) ^
               (((p_src[6 ^ ofs_h] >> (56 ^ ofs_l)) & 127) << 56);
            p_dest[7] =  ps_h[7] ^ ps_l[7] ^
               (((p_src[7 ^ ofs_h] >> (0 ^ ofs_l)) & 127) << 0) ^
               (((p_src[7 ^ ofs_h] >> (8 ^ ofs_l)) & 127) << 8) ^
               (((p_src[7 ^ ofs_h] >> (16 ^ ofs_l)) & 127) << 16) ^
               (((p_src[7 ^ ofs_h] >> (24 ^ ofs_l)) & 127) << 24) ^
               (((p_src[7 ^ ofs_h] >> (32 ^ ofs_l)) & 127) << 32) ^
               (((p_src[7 ^ ofs_h] >> (40 ^ ofs_l)) & 127) << 40) ^
               (((p_src[7 ^ ofs_h] >> (48 ^ ofs_l)) & 127) << 48) ^
               (((p_src[7 ^ ofs_h] >> (56 ^ ofs_l)) & 127) << 56);
            // %%END FOR
            p_src += 8; 
            p_dest += 8; 
            ++p_T;
        }
    }


    // Step 3: do rows with 24 entries, tags A, B, C // TODO: comment properly!!!!
    {
        uint_mmv_t mask[16];
        uint_mmv_t neg_mask[4];
        uint_mmv_t f = p_op->f_i, ef = p_op->ef_i, eps;
        for (i = 0; i < 3; ++i) {
             mask[i] = f >> (i << 3);
             mask[i + 4] = ef >> (i << 3);
        }
        neg_mask[0] = 0x7f7f7f7f7f7f7f7fULL;
        neg_mask[1] = 0x7f7f7f7f7f7f7f7fULL;
        neg_mask[2] = 0x7f7f7f7f7f7f7f7fULL;
        neg_mask[3] = 0;
        for (i = 0; i < 8; ++i) {
            uint_mmv_t x = mask[i];
            // %%MMV_UINT_SPREAD x, x
            // Spread bits 0,...,7 of x to the (8-bit long) fields
            // of x. A field of x is set to 0x7f if its 
            // corresponding bit in input x is one and to 0 otherwise.
            x = (x & 0xfULL) + ((x & 0xf0ULL) << 28);
            x = (x & 0x300000003ULL) 
                +  ((x & 0xc0000000cULL) << 14);
            x = (x & 0x1000100010001ULL) 
                +  ((x & 0x2000200020002ULL) << 7);
            x = (((x) << 7) - (x));
            // Bit spreading done.
            mask[i] = x = x & neg_mask[i & 3];
            mask[i + 8] = x ^ neg_mask[i & 3];
        }

        f <<= 3;
        ef <<= 3;
        eps = -(p_op->eps & 0x1ULL);
        for (i = 0; i < 96; i += 4) {
            uint_mmv_t t, t1, t2;
            // %%FOR j in range(V24_INTS_USED)
            // process uint_mmv_t 0 of row i/4 for tags A, B, C
            t1 = v_in[i + MM_OP127_OFS_A + 0]; 
            t = mask[0 + (f & 0x8ULL)];
            v_out[i + MM_OP127_OFS_A + 0] = t1 ^ t; 
            t1 = v_in[i + MM_OP127_OFS_B + 0]; 
            t2 = v_in[i + MM_OP127_OFS_C + 0];
            t &= (t1 ^ t2);
            t ^= mask[4 + (ef & 0x8ULL)];
            v_out[i + MM_OP127_OFS_B + 0] = t1 ^ t;
            t2 ^= eps & neg_mask[0];
            v_out[i + MM_OP127_OFS_C + 0] = t2 ^ t;
            // process uint_mmv_t 1 of row i/4 for tags A, B, C
            t1 = v_in[i + MM_OP127_OFS_A + 1]; 
            t = mask[1 + (f & 0x8ULL)];
            v_out[i + MM_OP127_OFS_A + 1] = t1 ^ t; 
            t1 = v_in[i + MM_OP127_OFS_B + 1]; 
            t2 = v_in[i + MM_OP127_OFS_C + 1];
            t &= (t1 ^ t2);
            t ^= mask[5 + (ef & 0x8ULL)];
            v_out[i + MM_OP127_OFS_B + 1] = t1 ^ t;
            t2 ^= eps & neg_mask[1];
            v_out[i + MM_OP127_OFS_C + 1] = t2 ^ t;
            // process uint_mmv_t 2 of row i/4 for tags A, B, C
            t1 = v_in[i + MM_OP127_OFS_A + 2]; 
            t = mask[2 + (f & 0x8ULL)];
            v_out[i + MM_OP127_OFS_A + 2] = t1 ^ t; 
            t1 = v_in[i + MM_OP127_OFS_B + 2]; 
            t2 = v_in[i + MM_OP127_OFS_C + 2];
            t &= (t1 ^ t2);
            t ^= mask[6 + (ef & 0x8ULL)];
            v_out[i + MM_OP127_OFS_B + 2] = t1 ^ t;
            t2 ^= eps & neg_mask[2];
            v_out[i + MM_OP127_OFS_C + 2] = t2 ^ t;
            // %%END FOR
            v_out[i + MM_OP127_OFS_A + 3] = 0;
            v_out[i + MM_OP127_OFS_B + 3] = 0;
            v_out[i + MM_OP127_OFS_C + 3] = 0;
            f >>= 1; ef >>= 1;      
        }
        //yet to be done!!!!
    }



    // If eps is odd: 
    //    negate entries X_d,i with scalar product <d,i> = 1
    if (p_op->eps & 1) mm127_neg_scalprod_d_i(v_out + MM_OP127_OFS_X); 
} 





// %%EXPORT p
void mm_op127_xy(uint_mmv_t *v_in, uint32_t f, uint32_t e, uint32_t eps, uint_mmv_t *v_out)
{
    mm_sub_op_xy_type s_op;
    mm_sub_prep_xy(f, e, eps, &s_op);
    mm_op127_do_xy(v_in, &s_op, v_out);
}


