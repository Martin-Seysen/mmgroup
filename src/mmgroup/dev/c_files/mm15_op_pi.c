/////////////////////////////////////////////////////////////////////////////
// This C file has been created automatically. Do not edit!!!
/////////////////////////////////////////////////////////////////////////////


// %%COMMENT
// TODO: Adjust this to new order of basis vectors rep!!!

#include <string.h>
#include "mm_op15.h"   
   
static const uint_mmv_t MM15_TBL_SCALPROD_HIGH[] = { 
// %%TABLE MMV_TBL_SCALPROD_HIGH, uint{INT_BITS}
0x0000000000000000ULL,0x0000000000000000ULL,
0xf0f00000f0f00000ULL,0x00000000ff000ff0ULL,
0x0ff000000ff00000ULL,0x00000000f0f0ff00ULL,
0xff000000ff000000ULL,0x000000000ff0f0f0ULL,
0xf0f000000000f0f0ULL,0x000000000ff0ff00ULL,
0x00000000f0f0f0f0ULL,0x00000000f0f0f0f0ULL,
0xff0000000ff0f0f0ULL,0x00000000ff000000ULL,
0x0ff00000ff00f0f0ULL,0x0000000000000ff0ULL,
0x0ff0000000000ff0ULL,0x00000000ff00f0f0ULL,
0xff000000f0f00ff0ULL,0x000000000000ff00ULL,
0x000000000ff00ff0ULL,0x000000000ff00ff0ULL,
0xf0f00000ff000ff0ULL,0x00000000f0f00000ULL,
0xff0000000000ff00ULL,0x00000000f0f00ff0ULL,
0x0ff00000f0f0ff00ULL,0x000000000ff00000ULL,
0xf0f000000ff0ff00ULL,0x000000000000f0f0ULL,
0x00000000ff00ff00ULL,0x00000000ff00ff00ULL,
0x000f000f000ffff0ULL,0x00000000000f000fULL,
0xf0ff000ff0fffff0ULL,0x00000000ff0f0fffULL,
0x0fff000f0ffffff0ULL,0x00000000f0ffff0fULL,
0xff0f000fff0ffff0ULL,0x000000000ffff0ffULL,
0xf0ff000f000f0f00ULL,0x000000000fffff0fULL,
0x000f000ff0ff0f00ULL,0x00000000f0fff0ffULL,
0xff0f000f0fff0f00ULL,0x00000000ff0f000fULL,
0x0fff000fff0f0f00ULL,0x00000000000f0fffULL,
0x0fff000f000ff000ULL,0x00000000ff0ff0ffULL,
0xff0f000ff0fff000ULL,0x00000000000fff0fULL,
0x000f000f0ffff000ULL,0x000000000fff0fffULL,
0xf0ff000fff0ff000ULL,0x00000000f0ff000fULL,
0xff0f000f000f00f0ULL,0x00000000f0ff0fffULL,
0x0fff000ff0ff00f0ULL,0x000000000fff000fULL,
0xf0ff000f0fff00f0ULL,0x00000000000ff0ffULL,
0x000f000fff0f00f0ULL,0x00000000ff0fff0fULL
};

static const uint_mmv_t MM15_TBL_SCALPROD_LOW[] = { 
// %%TABLE MMV_TBL_SCALPROD_LOW, uint{INT_BITS}
0x0000000000000000ULL,0x0000000000000000ULL,
0xffffffffffff0000ULL,0x00000000ffff0000ULL,
0xffffffffffff0000ULL,0x000000000000ffffULL,
0x0000000000000000ULL,0x00000000ffffffffULL,
0xf0f0f0f000000000ULL,0x00000000f0f0f0f0ULL,
0x0f0f0f0fffff0000ULL,0x000000000f0ff0f0ULL,
0x0f0f0f0fffff0000ULL,0x00000000f0f00f0fULL,
0xf0f0f0f000000000ULL,0x000000000f0f0f0fULL,
0x0ff00ff000000000ULL,0x000000000ff00ff0ULL,
0xf00ff00fffff0000ULL,0x00000000f00f0ff0ULL,
0xf00ff00fffff0000ULL,0x000000000ff0f00fULL,
0x0ff00ff000000000ULL,0x00000000f00ff00fULL,
0xff00ff0000000000ULL,0x00000000ff00ff00ULL,
0x00ff00ffffff0000ULL,0x0000000000ffff00ULL,
0x00ff00ffffff0000ULL,0x00000000ff0000ffULL,
0xff00ff0000000000ULL,0x0000000000ff00ffULL
};



// %%EXPORT
void mm15_neg_scalprod_d_i(uint_mmv_t* v)
// negate entries d (x) i with sclar product equal to 1
{
    const uint_mmv_t* p0 = MM15_TBL_SCALPROD_HIGH;
    const uint_mmv_t* p0_end = p0 + 32 * 2;

    // inversion of entries (d (x) i) with scalar product 1
    for (; p0 < p0_end; p0 += 2) {
        const uint_mmv_t* p1 = MM15_TBL_SCALPROD_LOW;
        const uint_mmv_t* p1_end = p1 + 16 * 2;
        for (; p1 < p1_end; p1 += 2) {
            // %%SCALAR_PROD_2048_UNROLL p0, p1, v 
            uint_mmv_t v_t;
            v[0] ^= (v_t = p0[0] ^ p1[0]);
            v[2] ^= v_t ^ 0xffff0000ffff0000ULL;
            v[4] ^= v_t ^ 0xffffffff0000ULL;
            v[6] ^= v_t ^ 0xffffffff00000000ULL;
            v[1] ^= (v_t = p0[1] ^ p1[1]);
            v[3] ^= v_t ^ 0xffffffffULL;
            v[5] ^= v_t ^ 0xffffffffULL;
            v[7] ^= v_t ^ 0x0ULL;
            v +=   4 * 2;
        }
    }
}





static uint32_t TABLE24_START[4] = {
   MM_OP15_OFS_X, MM_OP15_OFS_Z, MM_OP15_OFS_Y, MM_OP15_OFS_A
};



//TODO: Adjust this to new function mat24_op_all_autpl!!!

// %%EXPORT
void mm_op15_do_pi(uint_mmv_t *v_in, mm_sub_op_pi_type *p_op, uint_mmv_t *v_out)
{
    uint_fast32_t i;
    uint32_t table24_dest[4];

    // %%IF PERM24_USE_BENES_NET
    // The following mask is used by the actual permutation code
    // %%PERM24_BENES_DECLARE "v"
    uint_mmv_t v0, v1;
    uint_mmv_t benes_mask[15]; 

    // Prepare mask array from Benes network
    // %%PERM24_BENES_PREPARE "p_op->benes_net", benes_mask
    {
        uint_mmv_t tmp; 
        uint_fast8_t i;
        static uint8_t tbl[] = {
            // %%TABLE table_prepare_perm24, uint8_t
        0x00,0x01,0x02,0x03,0x10,0x11,0x12,0x04,
        0x05,0x06,0x07,0x08,0x16,0x17,0x18
        };

        for(i = 0; i < 15; ++i) {
            tmp = tbl[i]; tmp = (p_op->benes_net)[tmp & 15] >> (tmp & 0xf0);
            // %%MMV_UINT_SPREAD tmp, tmp
            // Spread bits 0,...,15 of tmp to the (4-bit long) fields
            // of tmp. A field of tmp is set to 0xf if its 
            // corresponding bit in input tmp is one and to 0 otherwise.
            tmp = (tmp & 0xffULL) + ((tmp & 0xff00ULL) << 24);
            tmp = (tmp & 0xf0000000fULL) 
                +  ((tmp & 0xf0000000f0ULL) << 12);
            tmp = (tmp & 0x3000300030003ULL) 
                +  ((tmp & 0xc000c000c000cULL) << 6);
            tmp = (tmp & 0x101010101010101ULL) 
                +  ((tmp & 0x202020202020202ULL) << 3);
            tmp = (((tmp) << 4) - (tmp));
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
            uint_mmv_t *ps = p_src + ((sgn_perm & 0x7ff) << 1);
            
            sgn_perm >>= (i + 12);  // sign for permutation
            // %%IF PERM24_USE_BENES_NET
            // Load 'ps' to temporary variables v0,...
            // %%PERM24_BENES_LOAD ps
            v0 = (ps)[0];
            v1 = (ps)[1];
            // Permute and possibly negate data in temp. variables
            // %%PERM24_BENES_PERMUTE benes_mask, sgn_perm
            // Permute the 24 small integers in '(v0, v1)' 
            // using the Benes network given by 'benes_mask'. All small   
            // integers are negated if bit 0 of 'sign' is set.
            sgn_perm = (-(sgn_perm & 0x1ULL)) & 0xffffffffffffffffULL;
            v0 ^= sgn_perm;
            sgn_perm &= 0xffffffffULL;
            v1 ^= sgn_perm;
            sgn_perm = (v0 ^ (v0 >> 4)) & benes_mask[0]; 
            v0 ^=  sgn_perm ^ (sgn_perm << 4);
            sgn_perm = (v0 ^ (v0 >> 8)) & benes_mask[1]; 
            v0 ^=  sgn_perm ^ (sgn_perm << 8);
            sgn_perm = (v0 ^ (v0 >> 16)) & benes_mask[2]; 
            v0 ^=  sgn_perm ^ (sgn_perm << 16);
            sgn_perm = (v0 ^ (v0 >> 32)) & benes_mask[3]; 
            v0 ^=  sgn_perm ^ (sgn_perm << 32);
            sgn_perm = (v1 ^ (v1 >> 4)) & benes_mask[4]; 
            v1 ^=  sgn_perm ^ (sgn_perm << 4);
            sgn_perm = (v1 ^ (v1 >> 8)) & benes_mask[5]; 
            v1 ^=  sgn_perm ^ (sgn_perm << 8);
            sgn_perm = (v1 ^ (v1 >> 16)) & benes_mask[6]; 
            v1 ^=  sgn_perm ^ (sgn_perm << 16);
            sgn_perm = (v0 ^ v1) & benes_mask[7]; 
            v0 ^=  sgn_perm;  v1 ^=  sgn_perm;
            sgn_perm = (v0 ^ (v0 >> 32)) & benes_mask[8]; 
            v0 ^=  sgn_perm ^ (sgn_perm << 32);
            sgn_perm = (v0 ^ (v0 >> 16)) & benes_mask[9]; 
            v0 ^=  sgn_perm ^ (sgn_perm << 16);
            sgn_perm = (v0 ^ (v0 >> 8)) & benes_mask[10]; 
            v0 ^=  sgn_perm ^ (sgn_perm << 8);
            sgn_perm = (v0 ^ (v0 >> 4)) & benes_mask[11]; 
            v0 ^=  sgn_perm ^ (sgn_perm << 4);
            sgn_perm = (v1 ^ (v1 >> 16)) & benes_mask[12]; 
            v1 ^=  sgn_perm ^ (sgn_perm << 16);
            sgn_perm = (v1 ^ (v1 >> 8)) & benes_mask[13]; 
            v1 ^=  sgn_perm ^ (sgn_perm << 8);
            sgn_perm = (v1 ^ (v1 >> 4)) & benes_mask[14]; 
            v1 ^=  sgn_perm ^ (sgn_perm << 4);
            // Permutation of small integers done.
            // Store temporary variables to 'p_dest'
            // %%PERM24_BENES_STORE p_dest
            (p_dest)[0] = v0 ;
            (p_dest)[1] = v1 ;
            // %%END IF
            p_dest +=  2;      
        }
        
    }    

    // Step 2: do rows with 64 entries // TODO: comment properly!!!!
    {
        // TODO: check this !!!!!!!!!!!!
        mm_sub_op_pi64_type *p_perm = p_op->tbl_perm64;
        uint8_t bytes[64];
        uint_mmv_t *p_out = v_out + MM_OP15_OFS_T;
        uint_mmv_t *p_end = p_out + 759 * 4;
        v_in +=  MM_OP15_OFS_T;
        for (; p_out < p_end; p_out += 4) {
            {
               uint_mmv_t v = p_perm->preimage;
               uint_mmv_t *p_in = v_in + ((v & 0x3ff) << 2);
               // %%LOAD_PERM64 p_in, bytes, v
               uint_mmv_t r0, r1;
               v = (-((v >> 12) & 1)) & 0xffffffffffffffffULL;
               r0 =  (p_in)[0] ^ (v);
               r1 = (r0 >> 4) & 0xf0f0f0f0f0f0f0fULL;
               r0 &= 0xf0f0f0f0f0f0f0fULL;
               (bytes)[0] = r0 >> 0;
               (bytes)[1] = r1 >> 0;
               (bytes)[2] = r0 >> 8;
               (bytes)[3] = r1 >> 8;
               (bytes)[4] = r0 >> 16;
               (bytes)[5] = r1 >> 16;
               (bytes)[6] = r0 >> 24;
               (bytes)[7] = r1 >> 24;
               (bytes)[8] = r0 >> 32;
               (bytes)[9] = r1 >> 32;
               (bytes)[10] = r0 >> 40;
               (bytes)[11] = r1 >> 40;
               (bytes)[12] = r0 >> 48;
               (bytes)[13] = r1 >> 48;
               (bytes)[14] = r0 >> 56;
               (bytes)[15] = r1 >> 56;
               r0 =  (p_in)[1] ^ (v);
               r1 = (r0 >> 4) & 0xf0f0f0f0f0f0f0fULL;
               r0 &= 0xf0f0f0f0f0f0f0fULL;
               (bytes)[16] = r0 >> 0;
               (bytes)[17] = r1 >> 0;
               (bytes)[18] = r0 >> 8;
               (bytes)[19] = r1 >> 8;
               (bytes)[20] = r0 >> 16;
               (bytes)[21] = r1 >> 16;
               (bytes)[22] = r0 >> 24;
               (bytes)[23] = r1 >> 24;
               (bytes)[24] = r0 >> 32;
               (bytes)[25] = r1 >> 32;
               (bytes)[26] = r0 >> 40;
               (bytes)[27] = r1 >> 40;
               (bytes)[28] = r0 >> 48;
               (bytes)[29] = r1 >> 48;
               (bytes)[30] = r0 >> 56;
               (bytes)[31] = r1 >> 56;
               r0 =  (p_in)[2] ^ (v);
               r1 = (r0 >> 4) & 0xf0f0f0f0f0f0f0fULL;
               r0 &= 0xf0f0f0f0f0f0f0fULL;
               (bytes)[32] = r0 >> 0;
               (bytes)[33] = r1 >> 0;
               (bytes)[34] = r0 >> 8;
               (bytes)[35] = r1 >> 8;
               (bytes)[36] = r0 >> 16;
               (bytes)[37] = r1 >> 16;
               (bytes)[38] = r0 >> 24;
               (bytes)[39] = r1 >> 24;
               (bytes)[40] = r0 >> 32;
               (bytes)[41] = r1 >> 32;
               (bytes)[42] = r0 >> 40;
               (bytes)[43] = r1 >> 40;
               (bytes)[44] = r0 >> 48;
               (bytes)[45] = r1 >> 48;
               (bytes)[46] = r0 >> 56;
               (bytes)[47] = r1 >> 56;
               r0 =  (p_in)[3] ^ (v);
               r1 = (r0 >> 4) & 0xf0f0f0f0f0f0f0fULL;
               r0 &= 0xf0f0f0f0f0f0f0fULL;
               (bytes)[48] = r0 >> 0;
               (bytes)[49] = r1 >> 0;
               (bytes)[50] = r0 >> 8;
               (bytes)[51] = r1 >> 8;
               (bytes)[52] = r0 >> 16;
               (bytes)[53] = r1 >> 16;
               (bytes)[54] = r0 >> 24;
               (bytes)[55] = r1 >> 24;
               (bytes)[56] = r0 >> 32;
               (bytes)[57] = r1 >> 32;
               (bytes)[58] = r0 >> 40;
               (bytes)[59] = r1 >> 40;
               (bytes)[60] = r0 >> 48;
               (bytes)[61] = r1 >> 48;
               (bytes)[62] = r0 >> 56;
               (bytes)[63] = r1 >> 56;
            }
            {
               // %%STORE_PERM64 bytes, p_out, "(p_perm->perm)"
               uint_mmv_t v;
               uint_fast8_t ri, r0, r1, r2, r3;
               r0 = (p_perm->perm)[0];
               r1 = (p_perm->perm)[1];
               r2 = (p_perm->perm)[2];
               r3 = (p_perm->perm)[3];
               ri = r0;
               v = (uint_mmv_t)(bytes[0]) << 0;
               v += (uint_mmv_t)(bytes[ri]) << 4;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 8;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 12;
               ri ^= r2;
               v += (uint_mmv_t)(bytes[ri]) << 16;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 20;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 24;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 28;
               ri ^= r3;
               v += (uint_mmv_t)(bytes[ri]) << 32;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 36;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 40;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 44;
               ri ^= r2;
               v += (uint_mmv_t)(bytes[ri]) << 48;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 52;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 56;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 60;
               p_out[0] = v ;
               ri ^= (p_perm->perm)[4];
               v = (uint_mmv_t)(bytes[ri]) << 0;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 4;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 8;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 12;
               ri ^= r2;
               v += (uint_mmv_t)(bytes[ri]) << 16;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 20;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 24;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 28;
               ri ^= r3;
               v += (uint_mmv_t)(bytes[ri]) << 32;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 36;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 40;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 44;
               ri ^= r2;
               v += (uint_mmv_t)(bytes[ri]) << 48;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 52;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 56;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 60;
               p_out[1] = v ;
               ri ^= (p_perm->perm)[5];
               v = (uint_mmv_t)(bytes[ri]) << 0;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 4;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 8;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 12;
               ri ^= r2;
               v += (uint_mmv_t)(bytes[ri]) << 16;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 20;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 24;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 28;
               ri ^= r3;
               v += (uint_mmv_t)(bytes[ri]) << 32;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 36;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 40;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 44;
               ri ^= r2;
               v += (uint_mmv_t)(bytes[ri]) << 48;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 52;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 56;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 60;
               p_out[2] = v ;
               ri ^= (p_perm->perm)[4];
               v = (uint_mmv_t)(bytes[ri]) << 0;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 4;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 8;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 12;
               ri ^= r2;
               v += (uint_mmv_t)(bytes[ri]) << 16;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 20;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 24;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 28;
               ri ^= r3;
               v += (uint_mmv_t)(bytes[ri]) << 32;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 36;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 40;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 44;
               ri ^= r2;
               v += (uint_mmv_t)(bytes[ri]) << 48;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 52;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 56;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 60;
               p_out[3] = v ;
            }
            ++p_perm;
        } 
    }

    // If d is odd: negate some entries    
    if (p_op->d & 1) {
        uint_fast16_t i;
        uint_mmv_t *v = v_out + MM_OP15_OFS_T;

        // Step odd 1:  negate suboctads of weight 4n+2 for tag T
        for (i = 0; i < 759; ++i) {
            // %%INVERT_PERM64 v
            v[0] ^= 0xf0fff0ffffff0ULL;
            v[1] ^= 0xf000000f000f0fffULL;
            v[2] ^= 0xf000000f000f0fffULL;
            v[3] ^= 0xfff0f000f000000fULL;
            v += 4;
        }

        mm15_neg_scalprod_d_i(v); 
    }
} 





// %%EXPORT p
void mm_op15_pi(uint_mmv_t *v_in, uint32_t delta, uint32_t pi, uint_mmv_t *v_out)
{
    mm_sub_op_pi_type s_op;
    mm_sub_prep_pi(delta, pi, &s_op);
    mm_op15_do_pi(v_in, &s_op, v_out);
}


// %%EXPORT p
void mm_op15_delta(uint_mmv_t *v_in, uint32_t delta, uint_mmv_t *v_out)
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
            sgn &= 0xffffffffffffffffULL;
            p_dest[0] = p_src[0]  ^ sgn;
            sgn &= 0xffffffffULL;
            p_dest[1] = p_src[1]  ^ sgn;
            // %%END FOR
            p_src +=  2;      
            p_dest +=  2;      
        }        
    }    

    // Step 2: do rows with 64 entries 
    // TODO: comment properly!!!!
    {
        v_in +=  MM_OP15_OFS_T;
        v_out += MM_OP15_OFS_T;
        for (i = 0; i < 759; ++i) {
            uint_mmv_t  sign = mat24_def_octad_to_gcode(i) & delta;
            sign ^=  sign >> 6; sign ^=  sign >> 3;
            sign = -((0x96 >> (sign & 7)) & 1);
            sign &= 0xffffffffffffffffULL;
            // %%FOR i in range({V64_INTS})
            v_out[0] = v_in[0]  ^  sign;
            v_out[1] = v_in[1]  ^  sign;
            v_out[2] = v_in[2]  ^  sign;
            v_out[3] = v_in[3]  ^  sign;
            // %%END FOR
            v_in += 4;
            v_out += 4;
        } 
        v_out -= 759 * 4 +  MM_OP15_OFS_T; // restore v_out
    }

    // If d is odd: negate some entries    
    if (delta & 1) {
        uint_fast16_t i;
        uint_mmv_t *v = v_out + MM_OP15_OFS_T;

        // Step odd 1:  negate suboctads of weight 4n+2 for tag T
        for (i = 0; i < 759; ++i) {
            // %%INVERT_PERM64 v
            v[0] ^= 0xf0fff0ffffff0ULL;
            v[1] ^= 0xf000000f000f0fffULL;
            v[2] ^= 0xf000000f000f0fffULL;
            v[3] ^= 0xfff0f000f000000fULL;
            v += 4;
        }

        mm15_neg_scalprod_d_i(v); 
    }
}


