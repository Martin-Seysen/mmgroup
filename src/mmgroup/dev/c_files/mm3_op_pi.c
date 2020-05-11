/////////////////////////////////////////////////////////////////////////////
// This C file has been created automatically. Do not edit!!!
/////////////////////////////////////////////////////////////////////////////


// %%COMMENT
// TODO: Adjust this to new order of basis vectors rep!!!

#include <string.h>
#include "mm_op3.h"   
   
static const uint_mmv_t MM3_TBL_SCALPROD_HIGH[] = { 
// %%TABLE MMV_TBL_SCALPROD_HIGH, uint{INT_BITS}
0x0000000000000000ULL,0x0000f03ccc00cc00ULL,
0x0000ccf03c003c00ULL,0x00003cccf000f000ULL,
0x00003cf0cc0000ccULL,0x0000cccc0000ccccULL,
0x0000f000f0003cccULL,0x0000003c3c00f0ccULL,
0x0000f0cc3c00003cULL,0x000000f0f000cc3cULL,
0x00003c3c00003c3cULL,0x0000cc00cc00f03cULL,
0x0000cc3cf00000f0ULL,0x00003c003c00ccf0ULL,
0x000000cccc003cf0ULL,0x0000f0f00000f0f0ULL,
0x00000303030303fcULL,0x0000f33fcf03cffcULL,
0x0000cff33f033ffcULL,0x00003fcff303f3fcULL,
0x00003ff3cf030330ULL,0x0000cfcf0303cf30ULL,
0x0000f303f3033f30ULL,0x0000033f3f03f330ULL,
0x0000f3cf3f0303c0ULL,0x000003f3f303cfc0ULL,
0x00003f3f03033fc0ULL,0x0000cf03cf03f3c0ULL,
0x0000cf3ff303030cULL,0x00003f033f03cf0cULL,
0x000003cfcf033f0cULL,0x0000f3f30303f30cULL
};

static const uint_mmv_t MM3_TBL_SCALPROD_LOW[] = { 
// %%TABLE MMV_TBL_SCALPROD_LOW, uint{INT_BITS}
0x0000000000000000ULL,0x0000ff00ffffff00ULL,
0x000000ffffffff00ULL,0x0000ffff00000000ULL,
0x0000cccccccc0000ULL,0x000033cc3333ff00ULL,
0x0000cc333333ff00ULL,0x00003333cccc0000ULL,
0x00003c3c3c3c0000ULL,0x0000c33cc3c3ff00ULL,
0x00003cc3c3c3ff00ULL,0x0000c3c33c3c0000ULL,
0x0000f0f0f0f00000ULL,0x00000ff00f0fff00ULL,
0x0000f00f0f0fff00ULL,0x00000f0ff0f00000ULL
};



// %%EXPORT
void mm3_neg_scalprod_d_i(uint_mmv_t* v)
// negate entries d (x) i with sclar product equal to 1
{
    const uint_mmv_t* p0 = MM3_TBL_SCALPROD_HIGH;
    const uint_mmv_t* p0_end = p0 + 32 * 1;

    // inversion of entries (d (x) i) with scalar product 1
    for (; p0 < p0_end; p0 += 1) {
        const uint_mmv_t* p1 = MM3_TBL_SCALPROD_LOW;
        const uint_mmv_t* p1_end = p1 + 16 * 1;
        for (; p1 < p1_end; p1 += 1) {
            // %%SCALAR_PROD_2048_UNROLL p0, p1, v 
            uint_mmv_t v_t;
            v[0] ^= (v_t = p0[0] ^ p1[0]);
            v[1] ^= v_t ^ 0xffffff00ff00ULL;
            v[2] ^= v_t ^ 0xffff00ffff00ULL;
            v[3] ^= v_t ^ 0xffff0000ULL;
            v +=   4 * 1;
        }
    }
}





static uint32_t TABLE24_START[4] = {
   MM_OP3_OFS_X, MM_OP3_OFS_Z, MM_OP3_OFS_Y, MM_OP3_OFS_A
};



//TODO: Adjust this to new function mat24_op_all_autpl!!!

// %%EXPORT
void mm_op3_do_pi(uint_mmv_t *v_in, mm_sub_op_pi_type *p_op, uint_mmv_t *v_out)
{
    uint_fast32_t i;
    uint32_t table24_dest[4];

    // %%IF PERM24_USE_BENES_NET
    // The following mask is used by the actual permutation code
    // %%PERM24_BENES_DECLARE "v"
    uint_mmv_t v0;
    uint_mmv_t benes_mask[9]; 

    // Prepare mask array from Benes network
    // %%PERM24_BENES_PREPARE "p_op->benes_net", benes_mask
    {
        uint_mmv_t tmp; 
        uint_fast8_t i;

        for(i = 0; i < 9; ++i) {
            tmp = (p_op->benes_net)[i];
            // %%MMV_UINT_SPREAD tmp, tmp
            // Spread bits 0,...,31 of tmp to the (2-bit long) fields
            // of tmp. A field of tmp is set to 0x3 if its 
            // corresponding bit in input tmp is one and to 0 otherwise.
            tmp = (tmp & 0xffffULL) 
                +  ((tmp & 0xffff0000ULL) << 16);
            tmp = (tmp & 0xff000000ffULL) 
                +  ((tmp & 0xff000000ff00ULL) << 8);
            tmp = (tmp & 0xf000f000f000fULL) 
                +  ((tmp & 0xf000f000f000f0ULL) << 4);
            tmp = (tmp & 0x303030303030303ULL) 
                +  ((tmp & 0xc0c0c0c0c0c0c0cULL) << 2);
            tmp = (tmp & 0x1111111111111111ULL) 
                +  ((tmp & 0x2222222222222222ULL) << 1);
            tmp = (((tmp) << 2) - (tmp));
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
            uint_mmv_t *ps = p_src + ((sgn_perm & 0x7ff) << 0);
            
            sgn_perm >>= (i + 12);  // sign for permutation
            // %%IF PERM24_USE_BENES_NET
            // Load 'ps' to temporary variables v0,...
            // %%PERM24_BENES_LOAD ps
            v0 = (ps)[0];
            // Permute and possibly negate data in temp. variables
            // %%PERM24_BENES_PERMUTE benes_mask, sgn_perm
            // Permute the 24 small integers in '(v0)' 
            // using the Benes network given by 'benes_mask'. All small   
            // integers are negated if bit 0 of 'sign' is set.
            sgn_perm = (-(sgn_perm & 0x1ULL)) & 0xffffffffffffULL;
            v0 ^= sgn_perm;
            sgn_perm = (v0 ^ (v0 >> 2)) & benes_mask[0]; 
            v0 ^=  sgn_perm ^ (sgn_perm << 2);
            sgn_perm = (v0 ^ (v0 >> 4)) & benes_mask[1]; 
            v0 ^=  sgn_perm ^ (sgn_perm << 4);
            sgn_perm = (v0 ^ (v0 >> 8)) & benes_mask[2]; 
            v0 ^=  sgn_perm ^ (sgn_perm << 8);
            sgn_perm = (v0 ^ (v0 >> 16)) & benes_mask[3]; 
            v0 ^=  sgn_perm ^ (sgn_perm << 16);
            sgn_perm = (v0 ^ (v0 >> 32)) & benes_mask[4]; 
            v0 ^=  sgn_perm ^ (sgn_perm << 32);
            sgn_perm = (v0 ^ (v0 >> 16)) & benes_mask[5]; 
            v0 ^=  sgn_perm ^ (sgn_perm << 16);
            sgn_perm = (v0 ^ (v0 >> 8)) & benes_mask[6]; 
            v0 ^=  sgn_perm ^ (sgn_perm << 8);
            sgn_perm = (v0 ^ (v0 >> 4)) & benes_mask[7]; 
            v0 ^=  sgn_perm ^ (sgn_perm << 4);
            sgn_perm = (v0 ^ (v0 >> 2)) & benes_mask[8]; 
            v0 ^=  sgn_perm ^ (sgn_perm << 2);
            // Permutation of small integers done.
            // Store temporary variables to 'p_dest'
            // %%PERM24_BENES_STORE p_dest
            (p_dest)[0] = v0 ;
            // %%END IF
            p_dest +=  1;      
        }
        
    }    

    // Step 2: do rows with 64 entries // TODO: comment properly!!!!
    {
        // TODO: check this !!!!!!!!!!!!
        mm_sub_op_pi64_type *p_perm = p_op->tbl_perm64;
        uint8_t bytes[64];
        uint_mmv_t *p_out = v_out + MM_OP3_OFS_T;
        uint_mmv_t *p_end = p_out + 759 * 2;
        v_in +=  MM_OP3_OFS_T;
        for (; p_out < p_end; p_out += 2) {
            {
               uint_mmv_t v = p_perm->preimage;
               uint_mmv_t *p_in = v_in + ((v & 0x3ff) << 1);
               // %%LOAD_PERM64 p_in, bytes, v
               uint_mmv_t r0, r1, r2, r3;
               v = (-((v >> 12) & 1)) & 0xffffffffffffffffULL;
               r0 =  (p_in)[0] ^ (v);
               r1 = (r0 >> 2) & 0x303030303030303ULL;
               r2 = (r0 >> 4) & 0x303030303030303ULL;
               r3 = (r0 >> 6) & 0x303030303030303ULL;
               r0 &= 0x303030303030303ULL;
               (bytes)[0] = r0 >> 0;
               (bytes)[1] = r1 >> 0;
               (bytes)[2] = r2 >> 0;
               (bytes)[3] = r3 >> 0;
               (bytes)[4] = r0 >> 8;
               (bytes)[5] = r1 >> 8;
               (bytes)[6] = r2 >> 8;
               (bytes)[7] = r3 >> 8;
               (bytes)[8] = r0 >> 16;
               (bytes)[9] = r1 >> 16;
               (bytes)[10] = r2 >> 16;
               (bytes)[11] = r3 >> 16;
               (bytes)[12] = r0 >> 24;
               (bytes)[13] = r1 >> 24;
               (bytes)[14] = r2 >> 24;
               (bytes)[15] = r3 >> 24;
               (bytes)[16] = r0 >> 32;
               (bytes)[17] = r1 >> 32;
               (bytes)[18] = r2 >> 32;
               (bytes)[19] = r3 >> 32;
               (bytes)[20] = r0 >> 40;
               (bytes)[21] = r1 >> 40;
               (bytes)[22] = r2 >> 40;
               (bytes)[23] = r3 >> 40;
               (bytes)[24] = r0 >> 48;
               (bytes)[25] = r1 >> 48;
               (bytes)[26] = r2 >> 48;
               (bytes)[27] = r3 >> 48;
               (bytes)[28] = r0 >> 56;
               (bytes)[29] = r1 >> 56;
               (bytes)[30] = r2 >> 56;
               (bytes)[31] = r3 >> 56;
               r0 =  (p_in)[1] ^ (v);
               r1 = (r0 >> 2) & 0x303030303030303ULL;
               r2 = (r0 >> 4) & 0x303030303030303ULL;
               r3 = (r0 >> 6) & 0x303030303030303ULL;
               r0 &= 0x303030303030303ULL;
               (bytes)[32] = r0 >> 0;
               (bytes)[33] = r1 >> 0;
               (bytes)[34] = r2 >> 0;
               (bytes)[35] = r3 >> 0;
               (bytes)[36] = r0 >> 8;
               (bytes)[37] = r1 >> 8;
               (bytes)[38] = r2 >> 8;
               (bytes)[39] = r3 >> 8;
               (bytes)[40] = r0 >> 16;
               (bytes)[41] = r1 >> 16;
               (bytes)[42] = r2 >> 16;
               (bytes)[43] = r3 >> 16;
               (bytes)[44] = r0 >> 24;
               (bytes)[45] = r1 >> 24;
               (bytes)[46] = r2 >> 24;
               (bytes)[47] = r3 >> 24;
               (bytes)[48] = r0 >> 32;
               (bytes)[49] = r1 >> 32;
               (bytes)[50] = r2 >> 32;
               (bytes)[51] = r3 >> 32;
               (bytes)[52] = r0 >> 40;
               (bytes)[53] = r1 >> 40;
               (bytes)[54] = r2 >> 40;
               (bytes)[55] = r3 >> 40;
               (bytes)[56] = r0 >> 48;
               (bytes)[57] = r1 >> 48;
               (bytes)[58] = r2 >> 48;
               (bytes)[59] = r3 >> 48;
               (bytes)[60] = r0 >> 56;
               (bytes)[61] = r1 >> 56;
               (bytes)[62] = r2 >> 56;
               (bytes)[63] = r3 >> 56;
            }
            {
               // %%STORE_PERM64 bytes, p_out, "(p_perm->perm)"
               uint_mmv_t v;
               uint_fast8_t ri, r0, r1, r2, r3, r4;
               r0 = (p_perm->perm)[0];
               r1 = (p_perm->perm)[1];
               r2 = (p_perm->perm)[2];
               r3 = (p_perm->perm)[3];
               r4 = (p_perm->perm)[4];
               ri = r0;
               v = (uint_mmv_t)(bytes[0]) << 0;
               v += (uint_mmv_t)(bytes[ri]) << 2;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 4;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 6;
               ri ^= r2;
               v += (uint_mmv_t)(bytes[ri]) << 8;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 10;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 12;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 14;
               ri ^= r3;
               v += (uint_mmv_t)(bytes[ri]) << 16;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 18;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 20;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 22;
               ri ^= r2;
               v += (uint_mmv_t)(bytes[ri]) << 24;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 26;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 28;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 30;
               ri ^= r4;
               v += (uint_mmv_t)(bytes[ri]) << 32;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 34;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 36;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 38;
               ri ^= r2;
               v += (uint_mmv_t)(bytes[ri]) << 40;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 42;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 44;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 46;
               ri ^= r3;
               v += (uint_mmv_t)(bytes[ri]) << 48;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 50;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 52;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 54;
               ri ^= r2;
               v += (uint_mmv_t)(bytes[ri]) << 56;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 58;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 60;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 62;
               p_out[0] = v ;
               ri ^= (p_perm->perm)[5];
               v = (uint_mmv_t)(bytes[ri]) << 0;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 2;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 4;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 6;
               ri ^= r2;
               v += (uint_mmv_t)(bytes[ri]) << 8;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 10;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 12;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 14;
               ri ^= r3;
               v += (uint_mmv_t)(bytes[ri]) << 16;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 18;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 20;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 22;
               ri ^= r2;
               v += (uint_mmv_t)(bytes[ri]) << 24;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 26;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 28;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 30;
               ri ^= r4;
               v += (uint_mmv_t)(bytes[ri]) << 32;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 34;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 36;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 38;
               ri ^= r2;
               v += (uint_mmv_t)(bytes[ri]) << 40;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 42;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 44;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 46;
               ri ^= r3;
               v += (uint_mmv_t)(bytes[ri]) << 48;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 50;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 52;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 54;
               ri ^= r2;
               v += (uint_mmv_t)(bytes[ri]) << 56;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 58;
               ri ^= r1;
               v += (uint_mmv_t)(bytes[ri]) << 60;
               ri ^= r0;
               v += (uint_mmv_t)(bytes[ri]) << 62;
               p_out[1] = v ;
            }
            ++p_perm;
        } 
    }

    // If d is odd: negate some entries    
    if (p_op->d & 1) {
        uint_fast16_t i;
        uint_mmv_t *v = v_out + MM_OP3_OFS_T;

        // Step odd 1:  negate suboctads of weight 4n+2 for tag T
        for (i = 0; i < 759; ++i) {
            // %%INVERT_PERM64 v
            v[0] ^= 0xc003033f033f3ffcULL;
            v[1] ^= 0xfcc0c003c003033fULL;
            v += 2;
        }

        mm3_neg_scalprod_d_i(v); 
    }
} 





// %%EXPORT p
void mm_op3_pi(uint_mmv_t *v_in, uint32_t delta, uint32_t pi, uint_mmv_t *v_out)
{
    mm_sub_op_pi_type s_op;
    mm_sub_prep_pi(delta, pi, &s_op);
    mm_op3_do_pi(v_in, &s_op, v_out);
}


// %%EXPORT p
void mm_op3_delta(uint_mmv_t *v_in, uint32_t delta, uint_mmv_t *v_out)
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
            // %%END FOR
            p_src +=  1;      
            p_dest +=  1;      
        }        
    }    

    // Step 2: do rows with 64 entries 
    // TODO: comment properly!!!!
    {
        v_in +=  MM_OP3_OFS_T;
        v_out += MM_OP3_OFS_T;
        for (i = 0; i < 759; ++i) {
            uint_mmv_t  sign = mat24_def_octad_to_gcode(i) & delta;
            sign ^=  sign >> 6; sign ^=  sign >> 3;
            sign = -((0x96 >> (sign & 7)) & 1);
            sign &= 0xffffffffffffffffULL;
            // %%FOR i in range({V64_INTS})
            v_out[0] = v_in[0]  ^  sign;
            v_out[1] = v_in[1]  ^  sign;
            // %%END FOR
            v_in += 2;
            v_out += 2;
        } 
        v_out -= 759 * 2 +  MM_OP3_OFS_T; // restore v_out
    }

    // If d is odd: negate some entries    
    if (delta & 1) {
        uint_fast16_t i;
        uint_mmv_t *v = v_out + MM_OP3_OFS_T;

        // Step odd 1:  negate suboctads of weight 4n+2 for tag T
        for (i = 0; i < 759; ++i) {
            // %%INVERT_PERM64 v
            v[0] ^= 0xc003033f033f3ffcULL;
            v[1] ^= 0xfcc0c003c003033fULL;
            v += 2;
        }

        mm3_neg_scalprod_d_i(v); 
    }
}


