/////////////////////////////////////////////////////////////////////////////
// This C file has been created automatically. Do not edit!!!
/////////////////////////////////////////////////////////////////////////////


// %%COMMENT
// TODO: Adjust this to new order of basis vectors rep!!!

#include <string.h>
#include "mm_op127.h"   
   
static const uint_mmv_t MM127_TBL_SCALPROD_HIGH[] = { 
// %%TABLE MMV_TBL_SCALPROD_HIGH, uint{INT_BITS}
0x0000000000000000ULL,0x0000000000000000ULL,
0x0000000000000000ULL,0x0000000000000000ULL,
0x7f007f0000000000ULL,0x7f007f0000000000ULL,
0x7f7f0000007f7f00ULL,0x0000000000000000ULL,
0x007f7f0000000000ULL,0x007f7f0000000000ULL,
0x7f007f007f7f0000ULL,0x0000000000000000ULL,
0x7f7f000000000000ULL,0x7f7f000000000000ULL,
0x007f7f007f007f00ULL,0x0000000000000000ULL,
0x000000007f007f00ULL,0x7f007f0000000000ULL,
0x007f7f007f7f0000ULL,0x0000000000000000ULL,
0x7f007f007f007f00ULL,0x0000000000000000ULL,
0x7f007f007f007f00ULL,0x0000000000000000ULL,
0x007f7f007f007f00ULL,0x7f7f000000000000ULL,
0x7f7f000000000000ULL,0x0000000000000000ULL,
0x7f7f00007f007f00ULL,0x007f7f0000000000ULL,
0x00000000007f7f00ULL,0x0000000000000000ULL,
0x00000000007f7f00ULL,0x007f7f0000000000ULL,
0x7f7f00007f007f00ULL,0x0000000000000000ULL,
0x7f007f00007f7f00ULL,0x7f7f000000000000ULL,
0x000000007f7f0000ULL,0x0000000000000000ULL,
0x007f7f00007f7f00ULL,0x0000000000000000ULL,
0x007f7f00007f7f00ULL,0x0000000000000000ULL,
0x7f7f0000007f7f00ULL,0x7f007f0000000000ULL,
0x7f007f0000000000ULL,0x0000000000000000ULL,
0x000000007f7f0000ULL,0x7f7f000000000000ULL,
0x7f007f00007f7f00ULL,0x0000000000000000ULL,
0x7f007f007f7f0000ULL,0x007f7f0000000000ULL,
0x007f7f0000000000ULL,0x0000000000000000ULL,
0x007f7f007f7f0000ULL,0x7f007f0000000000ULL,
0x000000007f007f00ULL,0x0000000000000000ULL,
0x7f7f00007f7f0000ULL,0x0000000000000000ULL,
0x7f7f00007f7f0000ULL,0x0000000000000000ULL,
0x0000007f7f7f7f00ULL,0x0000007f0000007fULL,
0x0000007f0000007fULL,0x0000000000000000ULL,
0x7f007f7f7f7f7f00ULL,0x7f007f7f0000007fULL,
0x7f7f007f007f7f7fULL,0x0000000000000000ULL,
0x007f7f7f7f7f7f00ULL,0x007f7f7f0000007fULL,
0x7f007f7f7f7f007fULL,0x0000000000000000ULL,
0x7f7f007f7f7f7f00ULL,0x7f7f007f0000007fULL,
0x007f7f7f7f007f7fULL,0x0000000000000000ULL,
0x0000007f007f0000ULL,0x7f007f7f0000007fULL,
0x007f7f7f7f7f007fULL,0x0000000000000000ULL,
0x7f007f7f007f0000ULL,0x0000007f0000007fULL,
0x7f007f7f7f007f7fULL,0x0000000000000000ULL,
0x007f7f7f007f0000ULL,0x7f7f007f0000007fULL,
0x7f7f007f0000007fULL,0x0000000000000000ULL,
0x7f7f007f007f0000ULL,0x007f7f7f0000007fULL,
0x0000007f007f7f7fULL,0x0000000000000000ULL,
0x0000007f7f000000ULL,0x007f7f7f0000007fULL,
0x7f7f007f7f007f7fULL,0x0000000000000000ULL,
0x7f007f7f7f000000ULL,0x7f7f007f0000007fULL,
0x0000007f7f7f007fULL,0x0000000000000000ULL,
0x007f7f7f7f000000ULL,0x0000007f0000007fULL,
0x007f7f7f007f7f7fULL,0x0000000000000000ULL,
0x7f7f007f7f000000ULL,0x7f007f7f0000007fULL,
0x7f007f7f0000007fULL,0x0000000000000000ULL,
0x0000007f00007f00ULL,0x7f7f007f0000007fULL,
0x7f007f7f007f7f7fULL,0x0000000000000000ULL,
0x7f007f7f00007f00ULL,0x007f7f7f0000007fULL,
0x007f7f7f0000007fULL,0x0000000000000000ULL,
0x007f7f7f00007f00ULL,0x7f007f7f0000007fULL,
0x0000007f7f007f7fULL,0x0000000000000000ULL,
0x7f7f007f00007f00ULL,0x0000007f0000007fULL,
0x7f7f007f7f7f007fULL,0x0000000000000000ULL
};

static const uint_mmv_t MM127_TBL_SCALPROD_LOW[] = { 
// %%TABLE MMV_TBL_SCALPROD_LOW, uint{INT_BITS}
0x0000000000000000ULL,0x0000000000000000ULL,
0x0000000000000000ULL,0x0000000000000000ULL,
0x7f7f7f7f00000000ULL,0x7f7f7f7f7f7f7f7fULL,
0x7f7f7f7f00000000ULL,0x0000000000000000ULL,
0x7f7f7f7f00000000ULL,0x7f7f7f7f7f7f7f7fULL,
0x000000007f7f7f7fULL,0x0000000000000000ULL,
0x0000000000000000ULL,0x0000000000000000ULL,
0x7f7f7f7f7f7f7f7fULL,0x0000000000000000ULL,
0x0000000000000000ULL,0x7f007f007f007f00ULL,
0x7f007f007f007f00ULL,0x0000000000000000ULL,
0x7f7f7f7f00000000ULL,0x007f007f007f007fULL,
0x007f007f7f007f00ULL,0x0000000000000000ULL,
0x7f7f7f7f00000000ULL,0x007f007f007f007fULL,
0x7f007f00007f007fULL,0x0000000000000000ULL,
0x0000000000000000ULL,0x7f007f007f007f00ULL,
0x007f007f007f007fULL,0x0000000000000000ULL,
0x0000000000000000ULL,0x007f7f00007f7f00ULL,
0x007f7f00007f7f00ULL,0x0000000000000000ULL,
0x7f7f7f7f00000000ULL,0x7f00007f7f00007fULL,
0x7f00007f007f7f00ULL,0x0000000000000000ULL,
0x7f7f7f7f00000000ULL,0x7f00007f7f00007fULL,
0x007f7f007f00007fULL,0x0000000000000000ULL,
0x0000000000000000ULL,0x007f7f00007f7f00ULL,
0x7f00007f7f00007fULL,0x0000000000000000ULL,
0x0000000000000000ULL,0x7f7f00007f7f0000ULL,
0x7f7f00007f7f0000ULL,0x0000000000000000ULL,
0x7f7f7f7f00000000ULL,0x00007f7f00007f7fULL,
0x00007f7f7f7f0000ULL,0x0000000000000000ULL,
0x7f7f7f7f00000000ULL,0x00007f7f00007f7fULL,
0x7f7f000000007f7fULL,0x0000000000000000ULL,
0x0000000000000000ULL,0x7f7f00007f7f0000ULL,
0x00007f7f00007f7fULL,0x0000000000000000ULL
};



// %%EXPORT
void mm127_neg_scalprod_d_i(uint_mmv_t* v)
// negate entries d (x) i with sclar product equal to 1
{
    const uint_mmv_t* p0 = MM127_TBL_SCALPROD_HIGH;
    const uint_mmv_t* p0_end = p0 + 32 * 4;

    // inversion of entries (d (x) i) with scalar product 1
    for (; p0 < p0_end; p0 += 4) {
        const uint_mmv_t* p1 = MM127_TBL_SCALPROD_LOW;
        const uint_mmv_t* p1_end = p1 + 16 * 4;
        for (; p1 < p1_end; p1 += 4) {
            // %%SCALAR_PROD_2048_UNROLL p0, p1, v 
            uint_mmv_t v_t;
            v[0] ^= (v_t = p0[0] ^ p1[0]);
            v[4] ^= v_t ^ 0x7f7f7f7f00000000ULL;
            v[8] ^= v_t ^ 0x7f7f7f7f00000000ULL;
            v[12] ^= v_t ^ 0x0ULL;
            v[1] ^= (v_t = p0[1] ^ p1[1]);
            v[5] ^= v_t ^ 0x7f7f7f7f00000000ULL;
            v[9] ^= v_t ^ 0x7f7f7f7fULL;
            v[13] ^= v_t ^ 0x7f7f7f7f7f7f7f7fULL;
            v[2] ^= (v_t = p0[2] ^ p1[2]);
            v[6] ^= v_t ^ 0x7f7f7f7f7f7f7f7fULL;
            v[10] ^= v_t ^ 0x7f7f7f7f7f7f7f7fULL;
            v[14] ^= v_t ^ 0x0ULL;
            v +=   4 * 4;
        }
    }
}





static uint32_t TABLE24_START[4] = {
   MM_OP127_OFS_X, MM_OP127_OFS_Z, MM_OP127_OFS_Y, MM_OP127_OFS_A
};



//TODO: Adjust this to new function mat24_op_all_autpl!!!

// %%EXPORT
void mm_op127_do_pi(uint_mmv_t *v_in, mm_sub_op_pi_type *p_op, uint_mmv_t *v_out)
{
    uint_fast32_t i;
    uint32_t table24_dest[4];

    // %%IF PERM24_USE_BENES_NET
    // The following mask is used by the actual permutation code
    // %%PERM24_BENES_DECLARE "v"
    uint_mmv_t v0, v1, v2;
    uint_mmv_t benes_mask[21]; 

    // Prepare mask array from Benes network
    // %%PERM24_BENES_PREPARE "p_op->benes_net", benes_mask
    {
        uint_mmv_t tmp; 
        uint_fast8_t i;
        static uint8_t tbl[] = {
            // %%TABLE table_prepare_perm24, uint8_t
        0x00,0x01,0x02,0x10,0x11,0x12,0x20,0x21,
        0x22,0x03,0x04,0x05,0x06,0x07,0x08,0x16,
        0x17,0x18,0x26,0x27,0x28
        };

        for(i = 0; i < 21; ++i) {
            tmp = tbl[i];
            tmp = (p_op->benes_net)[tmp & 15] >> ((tmp & 0xf0) >> 1);
            // %%MMV_UINT_SPREAD tmp, tmp
            // Spread bits 0,...,7 of tmp to the (8-bit long) fields
            // of tmp. A field of tmp is set to 0x7f if its 
            // corresponding bit in input tmp is one and to 0 otherwise.
            tmp = (tmp & 0xfULL) + ((tmp & 0xf0ULL) << 28);
            tmp = (tmp & 0x300000003ULL) 
                +  ((tmp & 0xc0000000cULL) << 14);
            tmp = (tmp & 0x1000100010001ULL) 
                +  ((tmp & 0x2000200020002ULL) << 7);
            tmp = (((tmp) << 7) - (tmp));
            // Bit spreading done.
            benes_mask[i] = tmp;
        }
    }
    // %%END IF
    
    // Step 1: do rows with 24 entries 
    // TODO: comment properly!!!!
    for (i = 0; i < 4; ++i) table24_dest[i] = TABLE24_START[i];
    i = (TABLE24_START[1] ^ TABLE24_START[2]) & -(p_op->d & 1);
    table24_dest[1] ^= i;  table24_dest[2] ^= i; 

    for (i = 0; i < 4; ++i)  {
        uint16_t *p_perm  = p_op->tbl_perm24_big +  (i < 3 ? 0 : 2048);
        uint_mmv_t *p_src = v_in + TABLE24_START[i];
        uint_mmv_t *p_dest = v_out + table24_dest[i];
        uint_fast32_t i1, len = i < 3 ? 2048 : 72;
        for (i1 = 0; i1 < len; ++i1) {
            uint_mmv_t sgn_perm = p_perm[i1];
            uint_mmv_t *ps = p_src + ((sgn_perm & 0x7ff) << 2);
            
            sgn_perm >>= (i + 12);  // sign for permutation
            // %%IF PERM24_USE_BENES_NET
            // Load 'ps' to temporary variables v0,...
            // %%PERM24_BENES_LOAD ps
            v0 = (ps)[0];
            v1 = (ps)[1];
            v2 = (ps)[2];
            // Permute and possibly negate data in temp. variables
            // %%PERM24_BENES_PERMUTE benes_mask, sgn_perm
            // Permute the 24 small integers in '(v0, v1, v2)' 
            // using the Benes network given by 'benes_mask'. All small   
            // integers are negated if bit 0 of 'sign' is set.
            sgn_perm = (-(sgn_perm & 0x1ULL)) & 0x7f7f7f7f7f7f7f7fULL;
            v0 ^= sgn_perm;
            v1 ^= sgn_perm;
            v2 ^= sgn_perm;
            sgn_perm = (v0 ^ (v0 >> 8)) & benes_mask[0]; 
            v0 ^=  sgn_perm ^ (sgn_perm << 8);
            sgn_perm = (v0 ^ (v0 >> 16)) & benes_mask[1]; 
            v0 ^=  sgn_perm ^ (sgn_perm << 16);
            sgn_perm = (v0 ^ (v0 >> 32)) & benes_mask[2]; 
            v0 ^=  sgn_perm ^ (sgn_perm << 32);
            sgn_perm = (v1 ^ (v1 >> 8)) & benes_mask[3]; 
            v1 ^=  sgn_perm ^ (sgn_perm << 8);
            sgn_perm = (v1 ^ (v1 >> 16)) & benes_mask[4]; 
            v1 ^=  sgn_perm ^ (sgn_perm << 16);
            sgn_perm = (v1 ^ (v1 >> 32)) & benes_mask[5]; 
            v1 ^=  sgn_perm ^ (sgn_perm << 32);
            sgn_perm = (v2 ^ (v2 >> 8)) & benes_mask[6]; 
            v2 ^=  sgn_perm ^ (sgn_perm << 8);
            sgn_perm = (v2 ^ (v2 >> 16)) & benes_mask[7]; 
            v2 ^=  sgn_perm ^ (sgn_perm << 16);
            sgn_perm = (v2 ^ (v2 >> 32)) & benes_mask[8]; 
            v2 ^=  sgn_perm ^ (sgn_perm << 32);
            sgn_perm = (v0 ^ v1) & benes_mask[9]; 
            v0 ^=  sgn_perm;  v1 ^=  sgn_perm;
            sgn_perm = (v0 ^ v2) & benes_mask[10]; 
            v0 ^=  sgn_perm;  v2 ^=  sgn_perm;
            sgn_perm = (v0 ^ v1) & benes_mask[11]; 
            v0 ^=  sgn_perm;  v1 ^=  sgn_perm;
            sgn_perm = (v0 ^ (v0 >> 32)) & benes_mask[12]; 
            v0 ^=  sgn_perm ^ (sgn_perm << 32);
            sgn_perm = (v0 ^ (v0 >> 16)) & benes_mask[13]; 
            v0 ^=  sgn_perm ^ (sgn_perm << 16);
            sgn_perm = (v0 ^ (v0 >> 8)) & benes_mask[14]; 
            v0 ^=  sgn_perm ^ (sgn_perm << 8);
            sgn_perm = (v1 ^ (v1 >> 32)) & benes_mask[15]; 
            v1 ^=  sgn_perm ^ (sgn_perm << 32);
            sgn_perm = (v1 ^ (v1 >> 16)) & benes_mask[16]; 
            v1 ^=  sgn_perm ^ (sgn_perm << 16);
            sgn_perm = (v1 ^ (v1 >> 8)) & benes_mask[17]; 
            v1 ^=  sgn_perm ^ (sgn_perm << 8);
            sgn_perm = (v2 ^ (v2 >> 32)) & benes_mask[18]; 
            v2 ^=  sgn_perm ^ (sgn_perm << 32);
            sgn_perm = (v2 ^ (v2 >> 16)) & benes_mask[19]; 
            v2 ^=  sgn_perm ^ (sgn_perm << 16);
            sgn_perm = (v2 ^ (v2 >> 8)) & benes_mask[20]; 
            v2 ^=  sgn_perm ^ (sgn_perm << 8);
            // Permutation of small integers done.
            // Store temporary variables to 'p_dest'
            // %%PERM24_BENES_STORE p_dest
            (p_dest)[0] = v0 ;
            (p_dest)[1] = v1 ;
            (p_dest)[2] = v2 ;
            // %%END IF
            p_dest +=  4;      
        }
        
    }    

    // Step 2: do rows with 64 entries // TODO: comment properly!!!!
    {
        // TODO: check this !!!!!!!!!!!!
        mm_sub_op_pi64_type *p_perm = p_op->tbl_perm64;
        uint8_t bytes[64];
        uint_mmv_t *p_out = v_out + MM_OP127_OFS_T;
        uint_mmv_t *p_end = p_out + 759 * 8;
        v_in +=  MM_OP127_OFS_T;
        for (; p_out < p_end; p_out += 8) {
            {
               uint_mmv_t v = p_perm->preimage;
               uint_mmv_t *p_in = v_in + ((v & 0x3ff) << 3);
               // %%LOAD_PERM64 p_in, bytes, v
               uint_mmv_t r0;
               v = (-((v >> 12) & 1)) & 0x7f7f7f7f7f7f7f7fULL;
               r0 =  (p_in)[0] ^ (v);
               r0 &= 0x7f7f7f7f7f7f7f7fULL;
               (bytes)[0] = r0 >> 0;
               (bytes)[1] = r0 >> 8;
               (bytes)[2] = r0 >> 16;
               (bytes)[3] = r0 >> 24;
               (bytes)[4] = r0 >> 32;
               (bytes)[5] = r0 >> 40;
               (bytes)[6] = r0 >> 48;
               (bytes)[7] = r0 >> 56;
               r0 =  (p_in)[1] ^ (v);
               r0 &= 0x7f7f7f7f7f7f7f7fULL;
               (bytes)[8] = r0 >> 0;
               (bytes)[9] = r0 >> 8;
               (bytes)[10] = r0 >> 16;
               (bytes)[11] = r0 >> 24;
               (bytes)[12] = r0 >> 32;
               (bytes)[13] = r0 >> 40;
               (bytes)[14] = r0 >> 48;
               (bytes)[15] = r0 >> 56;
               r0 =  (p_in)[2] ^ (v);
               r0 &= 0x7f7f7f7f7f7f7f7fULL;
               (bytes)[16] = r0 >> 0;
               (bytes)[17] = r0 >> 8;
               (bytes)[18] = r0 >> 16;
               (bytes)[19] = r0 >> 24;
               (bytes)[20] = r0 >> 32;
               (bytes)[21] = r0 >> 40;
               (bytes)[22] = r0 >> 48;
               (bytes)[23] = r0 >> 56;
               r0 =  (p_in)[3] ^ (v);
               r0 &= 0x7f7f7f7f7f7f7f7fULL;
               (bytes)[24] = r0 >> 0;
               (bytes)[25] = r0 >> 8;
               (bytes)[26] = r0 >> 16;
               (bytes)[27] = r0 >> 24;
               (bytes)[28] = r0 >> 32;
               (bytes)[29] = r0 >> 40;
               (bytes)[30] = r0 >> 48;
               (bytes)[31] = r0 >> 56;
               r0 =  (p_in)[4] ^ (v);
               r0 &= 0x7f7f7f7f7f7f7f7fULL;
               (bytes)[32] = r0 >> 0;
               (bytes)[33] = r0 >> 8;
               (bytes)[34] = r0 >> 16;
               (bytes)[35] = r0 >> 24;
               (bytes)[36] = r0 >> 32;
               (bytes)[37] = r0 >> 40;
               (bytes)[38] = r0 >> 48;
               (bytes)[39] = r0 >> 56;
               r0 =  (p_in)[5] ^ (v);
               r0 &= 0x7f7f7f7f7f7f7f7fULL;
               (bytes)[40] = r0 >> 0;
               (bytes)[41] = r0 >> 8;
               (bytes)[42] = r0 >> 16;
               (bytes)[43] = r0 >> 24;
               (bytes)[44] = r0 >> 32;
               (bytes)[45] = r0 >> 40;
               (bytes)[46] = r0 >> 48;
               (bytes)[47] = r0 >> 56;
               r0 =  (p_in)[6] ^ (v);
               r0 &= 0x7f7f7f7f7f7f7f7fULL;
               (bytes)[48] = r0 >> 0;
               (bytes)[49] = r0 >> 8;
               (bytes)[50] = r0 >> 16;
               (bytes)[51] = r0 >> 24;
               (bytes)[52] = r0 >> 32;
               (bytes)[53] = r0 >> 40;
               (bytes)[54] = r0 >> 48;
               (bytes)[55] = r0 >> 56;
               r0 =  (p_in)[7] ^ (v);
               r0 &= 0x7f7f7f7f7f7f7f7fULL;
               (bytes)[56] = r0 >> 0;
               (bytes)[57] = r0 >> 8;
               (bytes)[58] = r0 >> 16;
               (bytes)[59] = r0 >> 24;
               (bytes)[60] = r0 >> 32;
               (bytes)[61] = r0 >> 40;
               (bytes)[62] = r0 >> 48;
               (bytes)[63] = r0 >> 56;
            }
            {
               // %%STORE_PERM64 bytes, p_out, "(p_perm->perm)"
               uint_mmv_t v;
               uint_fast8_t ri, r0, r1, r2;
               r0 = (p_perm->perm)[0];
               r1 = (p_perm->perm)[1];
               r2 = (p_perm->perm)[2];
               ri = r0;
               v = (uint_mmv_t)(bytes[0]) << 0;
               v += (uint_mmv_t)(bytes[ri]) << 8;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 16;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 24;
               ri ^= r2;
               v += (uint_mmv_t)(bytes[ri]) << 32;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 40;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 48;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 56;
               p_out[0] = v ;
               ri ^= (p_perm->perm)[3];
               v = (uint_mmv_t)(bytes[ri]) << 0;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 8;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 16;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 24;
               ri ^= r2;
               v += (uint_mmv_t)(bytes[ri]) << 32;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 40;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 48;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 56;
               p_out[1] = v ;
               ri ^= (p_perm->perm)[4];
               v = (uint_mmv_t)(bytes[ri]) << 0;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 8;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 16;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 24;
               ri ^= r2;
               v += (uint_mmv_t)(bytes[ri]) << 32;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 40;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 48;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 56;
               p_out[2] = v ;
               ri ^= (p_perm->perm)[3];
               v = (uint_mmv_t)(bytes[ri]) << 0;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 8;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 16;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 24;
               ri ^= r2;
               v += (uint_mmv_t)(bytes[ri]) << 32;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 40;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 48;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 56;
               p_out[3] = v ;
               ri ^= (p_perm->perm)[5];
               v = (uint_mmv_t)(bytes[ri]) << 0;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 8;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 16;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 24;
               ri ^= r2;
               v += (uint_mmv_t)(bytes[ri]) << 32;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 40;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 48;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 56;
               p_out[4] = v ;
               ri ^= (p_perm->perm)[3];
               v = (uint_mmv_t)(bytes[ri]) << 0;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 8;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 16;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 24;
               ri ^= r2;
               v += (uint_mmv_t)(bytes[ri]) << 32;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 40;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 48;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 56;
               p_out[5] = v ;
               ri ^= (p_perm->perm)[4];
               v = (uint_mmv_t)(bytes[ri]) << 0;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 8;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 16;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 24;
               ri ^= r2;
               v += (uint_mmv_t)(bytes[ri]) << 32;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 40;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 48;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 56;
               p_out[6] = v ;
               ri ^= (p_perm->perm)[3];
               v = (uint_mmv_t)(bytes[ri]) << 0;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 8;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 16;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 24;
               ri ^= r2;
               v += (uint_mmv_t)(bytes[ri]) << 32;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 40;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 48;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 56;
               p_out[7] = v ;
            }
            ++p_perm;
        } 
    }

    // If d is odd: negate some entries    
    if (p_op->d & 1) {
        uint_fast16_t i;
        uint_mmv_t *v = v_out + MM_OP127_OFS_T;

        // Step odd 1:  negate suboctads of weight 4n+2 for tag T
        for (i = 0; i < 759; ++i) {
            // %%INVERT_PERM64 v
            v[0] ^= 0x7f7f7f7f7f7f00ULL;
            v[1] ^= 0x7f007f7f7fULL;
            v[2] ^= 0x7f007f7f7fULL;
            v[3] ^= 0x7f0000000000007fULL;
            v[4] ^= 0x7f007f7f7fULL;
            v[5] ^= 0x7f0000000000007fULL;
            v[6] ^= 0x7f0000000000007fULL;
            v[7] ^= 0x7f7f7f007f000000ULL;
            v += 8;
        }

        mm127_neg_scalprod_d_i(v); 
    }
} 





// %%EXPORT p
void mm_op127_pi(uint_mmv_t *v_in, uint32_t delta, uint32_t pi, uint_mmv_t *v_out)
{
    mm_sub_op_pi_type s_op;
    mm_sub_prep_pi(delta, pi, &s_op);
    mm_op127_do_pi(v_in, &s_op, v_out);
}


// %%EXPORT p
void mm_op127_delta(uint_mmv_t *v_in, uint32_t delta, uint_mmv_t *v_out)
{
    uint_fast32_t i;
    uint8_t signs[2048];
    uint32_t table24_dest[4];

    mat24_op_all_cocode(delta, signs);
    for (i = 0; i < 72; ++i) signs[i] &= 7;
    for (i = 48; i < 72; ++i) signs[i] |= delta << 3;
    for (i = 0; i < 4; ++i) table24_dest[i] = TABLE24_START[i];
    i = (TABLE24_START[1] ^ TABLE24_START[2]) & -(delta & 1);
    table24_dest[1] ^= i;  table24_dest[2] ^= i; 

    // Step 1: do rows with 24 entries 
    // TODO: comment properly!!!!
    for (i = 0; i < 4; ++i)  {
        uint_mmv_t *p_src = v_in + TABLE24_START[i];
        uint_mmv_t *p_dest = v_out + table24_dest[i];
        uint_fast32_t i1, len = i < 3 ? 2048 : 72;
        for (i1 = 0; i1 < len; ++i1) {
            uint_mmv_t sgn = -((signs[i1] >> i) & 1);
            // %%FOR i in range(V24_INTS_USED)
            sgn &= 0x7f7f7f7f7f7f7f7fULL;
            p_dest[0] = p_src[0]  ^ sgn;
            p_dest[1] = p_src[1]  ^ sgn;
            p_dest[2] = p_src[2]  ^ sgn;
            // %%END FOR
            p_dest[3] = 0;
            p_src +=  4;      
            p_dest +=  4;      
        }        
    }    

    // Step 2: do rows with 64 entries 
    // TODO: comment properly!!!!
    {
        v_in +=  MM_OP127_OFS_T;
        v_out += MM_OP127_OFS_T;
        for (i = 0; i < 759; ++i) {
            uint_mmv_t  sign = mat24_def_octad_to_gcode(i) & delta;
            sign ^=  sign >> 6; sign ^=  sign >> 3;
            sign = -((0x96 >> (sign & 7)) & 1);
            sign &= 0x7f7f7f7f7f7f7f7fULL;
            // %%FOR i in range({V64_INTS})
            v_out[0] = v_in[0]  ^  sign;
            v_out[1] = v_in[1]  ^  sign;
            v_out[2] = v_in[2]  ^  sign;
            v_out[3] = v_in[3]  ^  sign;
            v_out[4] = v_in[4]  ^  sign;
            v_out[5] = v_in[5]  ^  sign;
            v_out[6] = v_in[6]  ^  sign;
            v_out[7] = v_in[7]  ^  sign;
            // %%END FOR
            v_in += 8;
            v_out += 8;
        } 
        v_out -= 759 * 8 +  MM_OP127_OFS_T; // restore v_out
    }

    // If d is odd: negate some entries    
    if (delta & 1) {
        uint_fast16_t i;
        uint_mmv_t *v = v_out + MM_OP127_OFS_T;

        // Step odd 1:  negate suboctads of weight 4n+2 for tag T
        for (i = 0; i < 759; ++i) {
            // %%INVERT_PERM64 v
            v[0] ^= 0x7f7f7f7f7f7f00ULL;
            v[1] ^= 0x7f007f7f7fULL;
            v[2] ^= 0x7f007f7f7fULL;
            v[3] ^= 0x7f0000000000007fULL;
            v[4] ^= 0x7f007f7f7fULL;
            v[5] ^= 0x7f0000000000007fULL;
            v[6] ^= 0x7f0000000000007fULL;
            v[7] ^= 0x7f7f7f007f000000ULL;
            v += 8;
        }

        mm127_neg_scalprod_d_i(v); 
    }
}


