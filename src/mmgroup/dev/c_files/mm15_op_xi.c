/////////////////////////////////////////////////////////////////////////////
// This C file has been created automatically. Do not edit!!!
/////////////////////////////////////////////////////////////////////////////

// %%COMMENT
// TODO: Yet to be documented!!!


#include <string.h>
#include "mat24_functions.h"
#include "mm_op15.h"   


static void mm_op15_xi_mon(uint_mmv_t *v_in,  uint32_t exp1, uint_mmv_t *v_out)
{
    uint_mmv_t *p_src, *p_dest;
    uint_fast32_t i, j;
    uint_fast32_t diff = exp1 ? 2048 : 0;
    uint8_t b[2496], *p_b;
    mm_sub_table_xi_type *p_tables = mm_sub_table_xi[exp1];
    uint16_t *p_perm;
    uint32_t *p_sign;



    ///////////////////////////////////////////////////////////////
    // Map tag BC to tag BC.
    ///////////////////////////////////////////////////////////////
    p_src = v_in + 48;
    p_dest = v_out + 48;
    p_sign = p_tables[0].p_sign;
    p_perm = p_tables[0].p_perm;

    for (i = 0; i < 1; ++i) {
        p_b = b;
        for (j = 0; j < 78; ++j) {
           // %%OP_XI_LOAD p_src, p_b, {i_src.SHAPE[2]}
           uint_mmv_t r0, r1;
           r0 =  (p_src)[0];
           r1 = (r0 >> 4) & 0xf0f0f0f0f0f0f0fULL;
           r0 &= 0xf0f0f0f0f0f0f0fULL;
           (p_b)[0] = r0 >> 0;
           (p_b)[1] = r1 >> 0;
           (p_b)[2] = r0 >> 8;
           (p_b)[3] = r1 >> 8;
           (p_b)[4] = r0 >> 16;
           (p_b)[5] = r1 >> 16;
           (p_b)[6] = r0 >> 24;
           (p_b)[7] = r1 >> 24;
           (p_b)[8] = r0 >> 32;
           (p_b)[9] = r1 >> 32;
           (p_b)[10] = r0 >> 40;
           (p_b)[11] = r1 >> 40;
           (p_b)[12] = r0 >> 48;
           (p_b)[13] = r1 >> 48;
           (p_b)[14] = r0 >> 56;
           (p_b)[15] = r1 >> 56;
           r0 =  (p_src)[1];
           r1 = (r0 >> 4) & 0xf0f0f0f0f0f0f0fULL;
           r0 &= 0xf0f0f0f0f0f0f0fULL;
           (p_b)[16] = r0 >> 0;
           (p_b)[17] = r1 >> 0;
           (p_b)[18] = r0 >> 8;
           (p_b)[19] = r1 >> 8;
           (p_b)[20] = r0 >> 16;
           (p_b)[21] = r1 >> 16;
           (p_b)[22] = r0 >> 24;
           (p_b)[23] = r1 >> 24;
           (p_b)[24] = r0 >> 32;
           (p_b)[25] = r1 >> 32;
           (p_b)[26] = r0 >> 40;
           (p_b)[27] = r1 >> 40;
           (p_b)[28] = r0 >> 48;
           (p_b)[29] = r1 >> 48;
           (p_b)[30] = r0 >> 56;
           (p_b)[31] = r1 >> 56;
           p_src += 2;
           p_b += 32;
        }

        for (j = 0; j < 78; ++j) {
           // %%OP_XI_STORE b, p_perm, p_sign, p_dest, {i_dest.SHAPE[2]}
           uint_mmv_t r0, r1;
           r0 = ((uint_mmv_t)(b[p_perm[0]]) << 0)
             + ((uint_mmv_t)(b[p_perm[1]]) << 4)
             + ((uint_mmv_t)(b[p_perm[2]]) << 8)
             + ((uint_mmv_t)(b[p_perm[3]]) << 12)
             + ((uint_mmv_t)(b[p_perm[4]]) << 16)
             + ((uint_mmv_t)(b[p_perm[5]]) << 20)
             + ((uint_mmv_t)(b[p_perm[6]]) << 24)
             + ((uint_mmv_t)(b[p_perm[7]]) << 28)
             + ((uint_mmv_t)(b[p_perm[8]]) << 32)
             + ((uint_mmv_t)(b[p_perm[9]]) << 36)
             + ((uint_mmv_t)(b[p_perm[10]]) << 40)
             + ((uint_mmv_t)(b[p_perm[11]]) << 44)
             + ((uint_mmv_t)(b[p_perm[12]]) << 48)
             + ((uint_mmv_t)(b[p_perm[13]]) << 52)
             + ((uint_mmv_t)(b[p_perm[14]]) << 56)
             + ((uint_mmv_t)(b[p_perm[15]]) << 60);
           r1 = (p_sign)[0] >> 0;
           // Spread bits 0,...,15 of r1 to the (4-bit long) fields
           // of r1. A field of r1 is set to 0xf if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xffULL) + ((r1 & 0xff00ULL) << 24);
           r1 = (r1 & 0xf0000000fULL) 
               +  ((r1 & 0xf0000000f0ULL) << 12);
           r1 = (r1 & 0x3000300030003ULL) 
               +  ((r1 & 0xc000c000c000cULL) << 6);
           r1 = (r1 & 0x101010101010101ULL) 
               +  ((r1 & 0x202020202020202ULL) << 3);
           r1 = (((r1) << 4) - (r1));
           // Bit spreading done.
           (p_dest)[0] = r0 ^ r1;
           r0 = ((uint_mmv_t)(b[p_perm[16]]) << 0)
             + ((uint_mmv_t)(b[p_perm[17]]) << 4)
             + ((uint_mmv_t)(b[p_perm[18]]) << 8)
             + ((uint_mmv_t)(b[p_perm[19]]) << 12)
             + ((uint_mmv_t)(b[p_perm[20]]) << 16)
             + ((uint_mmv_t)(b[p_perm[21]]) << 20)
             + ((uint_mmv_t)(b[p_perm[22]]) << 24)
             + ((uint_mmv_t)(b[p_perm[23]]) << 28)
             + ((uint_mmv_t)(b[p_perm[24]]) << 32)
             + ((uint_mmv_t)(b[p_perm[25]]) << 36)
             + ((uint_mmv_t)(b[p_perm[26]]) << 40)
             + ((uint_mmv_t)(b[p_perm[27]]) << 44)
             + ((uint_mmv_t)(b[p_perm[28]]) << 48)
             + ((uint_mmv_t)(b[p_perm[29]]) << 52)
             + ((uint_mmv_t)(b[p_perm[30]]) << 56)
             + ((uint_mmv_t)(b[p_perm[31]]) << 60);
           r1 = (p_sign)[0] >> 16;
           // Spread bits 0,...,15 of r1 to the (4-bit long) fields
           // of r1. A field of r1 is set to 0xf if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xffULL) + ((r1 & 0xff00ULL) << 24);
           r1 = (r1 & 0xf0000000fULL) 
               +  ((r1 & 0xf0000000f0ULL) << 12);
           r1 = (r1 & 0x3000300030003ULL) 
               +  ((r1 & 0xc000c000c000cULL) << 6);
           r1 = (r1 & 0x101010101010101ULL) 
               +  ((r1 & 0x202020202020202ULL) << 3);
           r1 = (((r1) << 4) - (r1));
           // Bit spreading done.
           (p_dest)[1] = r0 ^ r1;
           p_dest += 2;
           p_perm += 32;
           p_sign += 1;
        }
    }

    ///////////////////////////////////////////////////////////////
    // Map tag T0 to tag T0.
    ///////////////////////////////////////////////////////////////
    p_src = v_in + 204;
    p_dest = v_out + 204;
    p_sign = p_tables[1].p_sign;
    p_perm = p_tables[1].p_perm;

    for (i = 0; i < 45; ++i) {
        p_b = b;
        for (j = 0; j < 16; ++j) {
           // %%OP_XI_LOAD p_src, p_b, {i_src.SHAPE[2]}
           uint_mmv_t r0, r1;
           r0 =  (p_src)[0];
           r1 = (r0 >> 4) & 0xf0f0f0f0f0f0f0fULL;
           r0 &= 0xf0f0f0f0f0f0f0fULL;
           (p_b)[0] = r0 >> 0;
           (p_b)[1] = r1 >> 0;
           (p_b)[2] = r0 >> 8;
           (p_b)[3] = r1 >> 8;
           (p_b)[4] = r0 >> 16;
           (p_b)[5] = r1 >> 16;
           (p_b)[6] = r0 >> 24;
           (p_b)[7] = r1 >> 24;
           (p_b)[8] = r0 >> 32;
           (p_b)[9] = r1 >> 32;
           (p_b)[10] = r0 >> 40;
           (p_b)[11] = r1 >> 40;
           (p_b)[12] = r0 >> 48;
           (p_b)[13] = r1 >> 48;
           (p_b)[14] = r0 >> 56;
           (p_b)[15] = r1 >> 56;
           r0 =  (p_src)[1];
           r1 = (r0 >> 4) & 0xf0f0f0f0f0f0f0fULL;
           r0 &= 0xf0f0f0f0f0f0f0fULL;
           (p_b)[16] = r0 >> 0;
           (p_b)[17] = r1 >> 0;
           (p_b)[18] = r0 >> 8;
           (p_b)[19] = r1 >> 8;
           (p_b)[20] = r0 >> 16;
           (p_b)[21] = r1 >> 16;
           (p_b)[22] = r0 >> 24;
           (p_b)[23] = r1 >> 24;
           (p_b)[24] = r0 >> 32;
           (p_b)[25] = r1 >> 32;
           (p_b)[26] = r0 >> 40;
           (p_b)[27] = r1 >> 40;
           (p_b)[28] = r0 >> 48;
           (p_b)[29] = r1 >> 48;
           (p_b)[30] = r0 >> 56;
           (p_b)[31] = r1 >> 56;
           p_src += 2;
           p_b += 32;
        }

        for (j = 0; j < 16; ++j) {
           // %%OP_XI_STORE b, p_perm, p_sign, p_dest, {i_dest.SHAPE[2]}
           uint_mmv_t r0, r1;
           r0 = ((uint_mmv_t)(b[p_perm[0]]) << 0)
             + ((uint_mmv_t)(b[p_perm[1]]) << 4)
             + ((uint_mmv_t)(b[p_perm[2]]) << 8)
             + ((uint_mmv_t)(b[p_perm[3]]) << 12)
             + ((uint_mmv_t)(b[p_perm[4]]) << 16)
             + ((uint_mmv_t)(b[p_perm[5]]) << 20)
             + ((uint_mmv_t)(b[p_perm[6]]) << 24)
             + ((uint_mmv_t)(b[p_perm[7]]) << 28)
             + ((uint_mmv_t)(b[p_perm[8]]) << 32)
             + ((uint_mmv_t)(b[p_perm[9]]) << 36)
             + ((uint_mmv_t)(b[p_perm[10]]) << 40)
             + ((uint_mmv_t)(b[p_perm[11]]) << 44)
             + ((uint_mmv_t)(b[p_perm[12]]) << 48)
             + ((uint_mmv_t)(b[p_perm[13]]) << 52)
             + ((uint_mmv_t)(b[p_perm[14]]) << 56)
             + ((uint_mmv_t)(b[p_perm[15]]) << 60);
           r1 = (p_sign)[0] >> 0;
           // Spread bits 0,...,15 of r1 to the (4-bit long) fields
           // of r1. A field of r1 is set to 0xf if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xffULL) + ((r1 & 0xff00ULL) << 24);
           r1 = (r1 & 0xf0000000fULL) 
               +  ((r1 & 0xf0000000f0ULL) << 12);
           r1 = (r1 & 0x3000300030003ULL) 
               +  ((r1 & 0xc000c000c000cULL) << 6);
           r1 = (r1 & 0x101010101010101ULL) 
               +  ((r1 & 0x202020202020202ULL) << 3);
           r1 = (((r1) << 4) - (r1));
           // Bit spreading done.
           (p_dest)[0] = r0 ^ r1;
           r0 = ((uint_mmv_t)(b[p_perm[16]]) << 0)
             + ((uint_mmv_t)(b[p_perm[17]]) << 4)
             + ((uint_mmv_t)(b[p_perm[18]]) << 8)
             + ((uint_mmv_t)(b[p_perm[19]]) << 12)
             + ((uint_mmv_t)(b[p_perm[20]]) << 16)
             + ((uint_mmv_t)(b[p_perm[21]]) << 20)
             + ((uint_mmv_t)(b[p_perm[22]]) << 24)
             + ((uint_mmv_t)(b[p_perm[23]]) << 28)
             + ((uint_mmv_t)(b[p_perm[24]]) << 32)
             + ((uint_mmv_t)(b[p_perm[25]]) << 36)
             + ((uint_mmv_t)(b[p_perm[26]]) << 40)
             + ((uint_mmv_t)(b[p_perm[27]]) << 44)
             + ((uint_mmv_t)(b[p_perm[28]]) << 48)
             + ((uint_mmv_t)(b[p_perm[29]]) << 52)
             + ((uint_mmv_t)(b[p_perm[30]]) << 56)
             + ((uint_mmv_t)(b[p_perm[31]]) << 60);
           r1 = (p_sign)[0] >> 16;
           // Spread bits 0,...,15 of r1 to the (4-bit long) fields
           // of r1. A field of r1 is set to 0xf if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xffULL) + ((r1 & 0xff00ULL) << 24);
           r1 = (r1 & 0xf0000000fULL) 
               +  ((r1 & 0xf0000000f0ULL) << 12);
           r1 = (r1 & 0x3000300030003ULL) 
               +  ((r1 & 0xc000c000c000cULL) << 6);
           r1 = (r1 & 0x101010101010101ULL) 
               +  ((r1 & 0x202020202020202ULL) << 3);
           r1 = (((r1) << 4) - (r1));
           // Bit spreading done.
           (p_dest)[1] = r0 ^ r1;
           p_dest += 2;
           p_perm += 32;
           p_sign += 1;
        }
    }

    ///////////////////////////////////////////////////////////////
    // Map tag T1 to tag X0 if e = 1
    // Map tag T1 to tag X1 if e = 2
    ///////////////////////////////////////////////////////////////
    p_src = v_in + 1644;
    p_dest = v_out + 3180;
    p_dest += diff;
    p_sign = p_tables[2].p_sign;
    p_perm = p_tables[2].p_perm;

    for (i = 0; i < 64; ++i) {
        p_b = b;
        for (j = 0; j < 12; ++j) {
           // %%OP_XI_LOAD p_src, p_b, {i_src.SHAPE[2]}
           uint_mmv_t r0, r1;
           r0 =  (p_src)[0];
           r1 = (r0 >> 4) & 0xf0f0f0f0f0f0f0fULL;
           r0 &= 0xf0f0f0f0f0f0f0fULL;
           (p_b)[0] = r0 >> 0;
           (p_b)[1] = r1 >> 0;
           (p_b)[2] = r0 >> 8;
           (p_b)[3] = r1 >> 8;
           (p_b)[4] = r0 >> 16;
           (p_b)[5] = r1 >> 16;
           (p_b)[6] = r0 >> 24;
           (p_b)[7] = r1 >> 24;
           (p_b)[8] = r0 >> 32;
           (p_b)[9] = r1 >> 32;
           (p_b)[10] = r0 >> 40;
           (p_b)[11] = r1 >> 40;
           (p_b)[12] = r0 >> 48;
           (p_b)[13] = r1 >> 48;
           (p_b)[14] = r0 >> 56;
           (p_b)[15] = r1 >> 56;
           r0 =  (p_src)[1];
           r1 = (r0 >> 4) & 0xf0f0f0f0f0f0f0fULL;
           r0 &= 0xf0f0f0f0f0f0f0fULL;
           (p_b)[16] = r0 >> 0;
           (p_b)[17] = r1 >> 0;
           (p_b)[18] = r0 >> 8;
           (p_b)[19] = r1 >> 8;
           (p_b)[20] = r0 >> 16;
           (p_b)[21] = r1 >> 16;
           (p_b)[22] = r0 >> 24;
           (p_b)[23] = r1 >> 24;
           (p_b)[24] = r0 >> 32;
           (p_b)[25] = r1 >> 32;
           (p_b)[26] = r0 >> 40;
           (p_b)[27] = r1 >> 40;
           (p_b)[28] = r0 >> 48;
           (p_b)[29] = r1 >> 48;
           (p_b)[30] = r0 >> 56;
           (p_b)[31] = r1 >> 56;
           p_src += 2;
           p_b += 32;
        }

        for (j = 0; j < 16; ++j) {
           // %%OP_XI_STORE b, p_perm, p_sign, p_dest, {i_dest.SHAPE[2]}
           uint_mmv_t r0, r1;
           r0 = ((uint_mmv_t)(b[p_perm[0]]) << 0)
             + ((uint_mmv_t)(b[p_perm[1]]) << 4)
             + ((uint_mmv_t)(b[p_perm[2]]) << 8)
             + ((uint_mmv_t)(b[p_perm[3]]) << 12)
             + ((uint_mmv_t)(b[p_perm[4]]) << 16)
             + ((uint_mmv_t)(b[p_perm[5]]) << 20)
             + ((uint_mmv_t)(b[p_perm[6]]) << 24)
             + ((uint_mmv_t)(b[p_perm[7]]) << 28)
             + ((uint_mmv_t)(b[p_perm[8]]) << 32)
             + ((uint_mmv_t)(b[p_perm[9]]) << 36)
             + ((uint_mmv_t)(b[p_perm[10]]) << 40)
             + ((uint_mmv_t)(b[p_perm[11]]) << 44)
             + ((uint_mmv_t)(b[p_perm[12]]) << 48)
             + ((uint_mmv_t)(b[p_perm[13]]) << 52)
             + ((uint_mmv_t)(b[p_perm[14]]) << 56)
             + ((uint_mmv_t)(b[p_perm[15]]) << 60);
           r1 = (p_sign)[0] >> 0;
           // Spread bits 0,...,15 of r1 to the (4-bit long) fields
           // of r1. A field of r1 is set to 0xf if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xffULL) + ((r1 & 0xff00ULL) << 24);
           r1 = (r1 & 0xf0000000fULL) 
               +  ((r1 & 0xf0000000f0ULL) << 12);
           r1 = (r1 & 0x3000300030003ULL) 
               +  ((r1 & 0xc000c000c000cULL) << 6);
           r1 = (r1 & 0x101010101010101ULL) 
               +  ((r1 & 0x202020202020202ULL) << 3);
           r1 = (((r1) << 4) - (r1));
           // Bit spreading done.
           (p_dest)[0] = r0 ^ r1;
           r0 = ((uint_mmv_t)(b[p_perm[16]]) << 0)
             + ((uint_mmv_t)(b[p_perm[17]]) << 4)
             + ((uint_mmv_t)(b[p_perm[18]]) << 8)
             + ((uint_mmv_t)(b[p_perm[19]]) << 12)
             + ((uint_mmv_t)(b[p_perm[20]]) << 16)
             + ((uint_mmv_t)(b[p_perm[21]]) << 20)
             + ((uint_mmv_t)(b[p_perm[22]]) << 24)
             + ((uint_mmv_t)(b[p_perm[23]]) << 28);
           r1 = (p_sign)[0] >> 16;
           // Spread bits 0,...,15 of r1 to the (4-bit long) fields
           // of r1. A field of r1 is set to 0xf if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xffULL) + ((r1 & 0xff00ULL) << 24);
           r1 = (r1 & 0xf0000000fULL) 
               +  ((r1 & 0xf0000000f0ULL) << 12);
           r1 = (r1 & 0x3000300030003ULL) 
               +  ((r1 & 0xc000c000c000cULL) << 6);
           r1 = (r1 & 0x101010101010101ULL) 
               +  ((r1 & 0x202020202020202ULL) << 3);
           r1 = (((r1) << 4) - (r1));
           // Bit spreading done.
           (p_dest)[1] = r0 ^ r1;
           p_dest += 2;
           p_perm += 24;
           p_sign += 1;
        }
    }

    ///////////////////////////////////////////////////////////////
    // Map tag X0 to tag X1 if e = 1
    // Map tag X1 to tag X0 if e = 2
    ///////////////////////////////////////////////////////////////
    p_src = v_in + 3180;
    p_src += diff;
    p_dest = v_out + 5228;
    p_dest -= diff;
    p_sign = p_tables[3].p_sign;
    p_perm = p_tables[3].p_perm;

    for (i = 0; i < 64; ++i) {
        p_b = b;
        for (j = 0; j < 16; ++j) {
           // %%OP_XI_LOAD p_src, p_b, {i_src.SHAPE[2]}
           uint_mmv_t r0, r1;
           r0 =  (p_src)[0];
           r1 = (r0 >> 4) & 0xf0f0f0f0f0f0f0fULL;
           r0 &= 0xf0f0f0f0f0f0f0fULL;
           (p_b)[0] = r0 >> 0;
           (p_b)[1] = r1 >> 0;
           (p_b)[2] = r0 >> 8;
           (p_b)[3] = r1 >> 8;
           (p_b)[4] = r0 >> 16;
           (p_b)[5] = r1 >> 16;
           (p_b)[6] = r0 >> 24;
           (p_b)[7] = r1 >> 24;
           (p_b)[8] = r0 >> 32;
           (p_b)[9] = r1 >> 32;
           (p_b)[10] = r0 >> 40;
           (p_b)[11] = r1 >> 40;
           (p_b)[12] = r0 >> 48;
           (p_b)[13] = r1 >> 48;
           (p_b)[14] = r0 >> 56;
           (p_b)[15] = r1 >> 56;
           r0 =  (p_src)[1];
           r1 = (r0 >> 4) & 0xf0f0f0f0f0f0f0fULL;
           r0 &= 0xf0f0f0f0f0f0f0fULL;
           (p_b)[16] = r0 >> 0;
           (p_b)[17] = r1 >> 0;
           (p_b)[18] = r0 >> 8;
           (p_b)[19] = r1 >> 8;
           (p_b)[20] = r0 >> 16;
           (p_b)[21] = r1 >> 16;
           (p_b)[22] = r0 >> 24;
           (p_b)[23] = r1 >> 24;
           p_src += 2;
           p_b += 32;
        }

        for (j = 0; j < 16; ++j) {
           // %%OP_XI_STORE b, p_perm, p_sign, p_dest, {i_dest.SHAPE[2]}
           uint_mmv_t r0, r1;
           r0 = ((uint_mmv_t)(b[p_perm[0]]) << 0)
             + ((uint_mmv_t)(b[p_perm[1]]) << 4)
             + ((uint_mmv_t)(b[p_perm[2]]) << 8)
             + ((uint_mmv_t)(b[p_perm[3]]) << 12)
             + ((uint_mmv_t)(b[p_perm[4]]) << 16)
             + ((uint_mmv_t)(b[p_perm[5]]) << 20)
             + ((uint_mmv_t)(b[p_perm[6]]) << 24)
             + ((uint_mmv_t)(b[p_perm[7]]) << 28)
             + ((uint_mmv_t)(b[p_perm[8]]) << 32)
             + ((uint_mmv_t)(b[p_perm[9]]) << 36)
             + ((uint_mmv_t)(b[p_perm[10]]) << 40)
             + ((uint_mmv_t)(b[p_perm[11]]) << 44)
             + ((uint_mmv_t)(b[p_perm[12]]) << 48)
             + ((uint_mmv_t)(b[p_perm[13]]) << 52)
             + ((uint_mmv_t)(b[p_perm[14]]) << 56)
             + ((uint_mmv_t)(b[p_perm[15]]) << 60);
           r1 = (p_sign)[0] >> 0;
           // Spread bits 0,...,15 of r1 to the (4-bit long) fields
           // of r1. A field of r1 is set to 0xf if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xffULL) + ((r1 & 0xff00ULL) << 24);
           r1 = (r1 & 0xf0000000fULL) 
               +  ((r1 & 0xf0000000f0ULL) << 12);
           r1 = (r1 & 0x3000300030003ULL) 
               +  ((r1 & 0xc000c000c000cULL) << 6);
           r1 = (r1 & 0x101010101010101ULL) 
               +  ((r1 & 0x202020202020202ULL) << 3);
           r1 = (((r1) << 4) - (r1));
           // Bit spreading done.
           (p_dest)[0] = r0 ^ r1;
           r0 = ((uint_mmv_t)(b[p_perm[16]]) << 0)
             + ((uint_mmv_t)(b[p_perm[17]]) << 4)
             + ((uint_mmv_t)(b[p_perm[18]]) << 8)
             + ((uint_mmv_t)(b[p_perm[19]]) << 12)
             + ((uint_mmv_t)(b[p_perm[20]]) << 16)
             + ((uint_mmv_t)(b[p_perm[21]]) << 20)
             + ((uint_mmv_t)(b[p_perm[22]]) << 24)
             + ((uint_mmv_t)(b[p_perm[23]]) << 28);
           r1 = (p_sign)[0] >> 16;
           // Spread bits 0,...,15 of r1 to the (4-bit long) fields
           // of r1. A field of r1 is set to 0xf if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xffULL) + ((r1 & 0xff00ULL) << 24);
           r1 = (r1 & 0xf0000000fULL) 
               +  ((r1 & 0xf0000000f0ULL) << 12);
           r1 = (r1 & 0x3000300030003ULL) 
               +  ((r1 & 0xc000c000c000cULL) << 6);
           r1 = (r1 & 0x101010101010101ULL) 
               +  ((r1 & 0x202020202020202ULL) << 3);
           r1 = (((r1) << 4) - (r1));
           // Bit spreading done.
           (p_dest)[1] = r0 ^ r1;
           p_dest += 2;
           p_perm += 24;
           p_sign += 1;
        }
    }

    ///////////////////////////////////////////////////////////////
    // Map tag X1 to tag T1 if e = 1
    // Map tag X0 to tag T1 if e = 2
    ///////////////////////////////////////////////////////////////
    p_src = v_in + 5228;
    p_src -= diff;
    p_dest = v_out + 1644;
    p_sign = p_tables[4].p_sign;
    p_perm = p_tables[4].p_perm;

    for (i = 0; i < 64; ++i) {
        p_b = b;
        for (j = 0; j < 16; ++j) {
           // %%OP_XI_LOAD p_src, p_b, {i_src.SHAPE[2]}
           uint_mmv_t r0, r1;
           r0 =  (p_src)[0];
           r1 = (r0 >> 4) & 0xf0f0f0f0f0f0f0fULL;
           r0 &= 0xf0f0f0f0f0f0f0fULL;
           (p_b)[0] = r0 >> 0;
           (p_b)[1] = r1 >> 0;
           (p_b)[2] = r0 >> 8;
           (p_b)[3] = r1 >> 8;
           (p_b)[4] = r0 >> 16;
           (p_b)[5] = r1 >> 16;
           (p_b)[6] = r0 >> 24;
           (p_b)[7] = r1 >> 24;
           (p_b)[8] = r0 >> 32;
           (p_b)[9] = r1 >> 32;
           (p_b)[10] = r0 >> 40;
           (p_b)[11] = r1 >> 40;
           (p_b)[12] = r0 >> 48;
           (p_b)[13] = r1 >> 48;
           (p_b)[14] = r0 >> 56;
           (p_b)[15] = r1 >> 56;
           r0 =  (p_src)[1];
           r1 = (r0 >> 4) & 0xf0f0f0f0f0f0f0fULL;
           r0 &= 0xf0f0f0f0f0f0f0fULL;
           (p_b)[16] = r0 >> 0;
           (p_b)[17] = r1 >> 0;
           (p_b)[18] = r0 >> 8;
           (p_b)[19] = r1 >> 8;
           (p_b)[20] = r0 >> 16;
           (p_b)[21] = r1 >> 16;
           (p_b)[22] = r0 >> 24;
           (p_b)[23] = r1 >> 24;
           p_src += 2;
           p_b += 32;
        }

        for (j = 0; j < 12; ++j) {
           // %%OP_XI_STORE b, p_perm, p_sign, p_dest, {i_dest.SHAPE[2]}
           uint_mmv_t r0, r1;
           r0 = ((uint_mmv_t)(b[p_perm[0]]) << 0)
             + ((uint_mmv_t)(b[p_perm[1]]) << 4)
             + ((uint_mmv_t)(b[p_perm[2]]) << 8)
             + ((uint_mmv_t)(b[p_perm[3]]) << 12)
             + ((uint_mmv_t)(b[p_perm[4]]) << 16)
             + ((uint_mmv_t)(b[p_perm[5]]) << 20)
             + ((uint_mmv_t)(b[p_perm[6]]) << 24)
             + ((uint_mmv_t)(b[p_perm[7]]) << 28)
             + ((uint_mmv_t)(b[p_perm[8]]) << 32)
             + ((uint_mmv_t)(b[p_perm[9]]) << 36)
             + ((uint_mmv_t)(b[p_perm[10]]) << 40)
             + ((uint_mmv_t)(b[p_perm[11]]) << 44)
             + ((uint_mmv_t)(b[p_perm[12]]) << 48)
             + ((uint_mmv_t)(b[p_perm[13]]) << 52)
             + ((uint_mmv_t)(b[p_perm[14]]) << 56)
             + ((uint_mmv_t)(b[p_perm[15]]) << 60);
           r1 = (p_sign)[0] >> 0;
           // Spread bits 0,...,15 of r1 to the (4-bit long) fields
           // of r1. A field of r1 is set to 0xf if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xffULL) + ((r1 & 0xff00ULL) << 24);
           r1 = (r1 & 0xf0000000fULL) 
               +  ((r1 & 0xf0000000f0ULL) << 12);
           r1 = (r1 & 0x3000300030003ULL) 
               +  ((r1 & 0xc000c000c000cULL) << 6);
           r1 = (r1 & 0x101010101010101ULL) 
               +  ((r1 & 0x202020202020202ULL) << 3);
           r1 = (((r1) << 4) - (r1));
           // Bit spreading done.
           (p_dest)[0] = r0 ^ r1;
           r0 = ((uint_mmv_t)(b[p_perm[16]]) << 0)
             + ((uint_mmv_t)(b[p_perm[17]]) << 4)
             + ((uint_mmv_t)(b[p_perm[18]]) << 8)
             + ((uint_mmv_t)(b[p_perm[19]]) << 12)
             + ((uint_mmv_t)(b[p_perm[20]]) << 16)
             + ((uint_mmv_t)(b[p_perm[21]]) << 20)
             + ((uint_mmv_t)(b[p_perm[22]]) << 24)
             + ((uint_mmv_t)(b[p_perm[23]]) << 28)
             + ((uint_mmv_t)(b[p_perm[24]]) << 32)
             + ((uint_mmv_t)(b[p_perm[25]]) << 36)
             + ((uint_mmv_t)(b[p_perm[26]]) << 40)
             + ((uint_mmv_t)(b[p_perm[27]]) << 44)
             + ((uint_mmv_t)(b[p_perm[28]]) << 48)
             + ((uint_mmv_t)(b[p_perm[29]]) << 52)
             + ((uint_mmv_t)(b[p_perm[30]]) << 56)
             + ((uint_mmv_t)(b[p_perm[31]]) << 60);
           r1 = (p_sign)[0] >> 16;
           // Spread bits 0,...,15 of r1 to the (4-bit long) fields
           // of r1. A field of r1 is set to 0xf if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xffULL) + ((r1 & 0xff00ULL) << 24);
           r1 = (r1 & 0xf0000000fULL) 
               +  ((r1 & 0xf0000000f0ULL) << 12);
           r1 = (r1 & 0x3000300030003ULL) 
               +  ((r1 & 0xc000c000c000cULL) << 6);
           r1 = (r1 & 0x101010101010101ULL) 
               +  ((r1 & 0x202020202020202ULL) << 3);
           r1 = (((r1) << 4) - (r1));
           // Bit spreading done.
           (p_dest)[1] = r0 ^ r1;
           p_dest += 2;
           p_perm += 32;
           p_sign += 1;
        }
    }
}

static uint_mmv_t TAB15_XI64_MASK[] = {
// %%TABLE TABLE_MUL_MATRIX_XI64, uint{INT_BITS}
0x000f000f000f000fULL,0x0000000000000000ULL,
0x000f000f000f000fULL,0xffffffffffffffffULL,
0x0000000000000000ULL,0x000f000f000f000fULL,
0xffffffffffffffffULL,0x000f000f000f000fULL
};


#define HALF_YZ_SHIFT 11

static uint32_t TAB15_XI64_OFFSET[2][4] = {
    {
        MM_OP15_OFS_Z + (2 << HALF_YZ_SHIFT),
        MM_OP15_OFS_Z + (0 << HALF_YZ_SHIFT),
        MM_OP15_OFS_Z + (1 << HALF_YZ_SHIFT),
        MM_OP15_OFS_Z + (3 << HALF_YZ_SHIFT)
    },
    {
        MM_OP15_OFS_Z + (1 << HALF_YZ_SHIFT),
        MM_OP15_OFS_Z + (2 << HALF_YZ_SHIFT),
        MM_OP15_OFS_Z + (0 << HALF_YZ_SHIFT),
        MM_OP15_OFS_Z + (3 << HALF_YZ_SHIFT)
    },
};




static void mm_op15_xi_yz(uint_mmv_t *v_in,  uint32_t exp1, uint_mmv_t *v_out)
{
    uint_fast32_t i;
    uint_mmv_t *p_mask =  TAB15_XI64_MASK + exp1;
    for (i = 0; i < 64; ++i) {
        // %%MUL_MATRIX_XI64 v_in, p_mask, v_out

        // This is an automatically generated matrix operation, do not change!
        {
        uint_mmv_t r0, r1, r2, r3, r4;
        uint_mmv_t r5, r6, r7, r8, r9;
        uint_mmv_t r10, r11, r12, r13, r14;
        uint_mmv_t r15, r16;

        uint_mmv_t a[48];
        uint_fast32_t i;

        // TODO: write comment!!!
        // 
        r0 = v_in[0] ^  p_mask[2];
        r16 = ((r0 ^ (r0 >> 4)) & 0xf000f000f000f0ULL);
        r0 ^= (r16 | (r16 << 4));
        a[16] = ((r0 >> 32) & 0xffffffffULL);
        r0 = (r0 & 0xffffffffULL);
        a[32] = v_in[1] ^  p_mask[2];
        r16 = ((a[32] ^ (a[32] >> 4))
            & 0xf000f000f000f0ULL);
        a[32] ^= (r16 | (r16 << 4));
        r1 = v_in[28] ^  p_mask[0];
        r16 = ((r1 ^ (r1 >> 4)) & 0xf000f000f000f0ULL);
        r1 ^= (r16 | (r16 << 4));
        a[17] = ((r1 >> 32) & 0xffffffffULL);
        r1 = (r1 & 0xffffffffULL);
        a[33] = v_in[29] ^  p_mask[0];
        r16 = ((a[33] ^ (a[33] >> 4))
            & 0xf000f000f000f0ULL);
        a[33] ^= (r16 | (r16 << 4));
        r2 = v_in[26] ^  p_mask[0];
        r16 = ((r2 ^ (r2 >> 4)) & 0xf000f000f000f0ULL);
        r2 ^= (r16 | (r16 << 4));
        a[18] = ((r2 >> 32) & 0xffffffffULL);
        r2 = (r2 & 0xffffffffULL);
        a[34] = v_in[27] ^  p_mask[0];
        r16 = ((a[34] ^ (a[34] >> 4))
            & 0xf000f000f000f0ULL);
        a[34] ^= (r16 | (r16 << 4));
        r3 = v_in[6] ^  p_mask[0];
        r16 = ((r3 ^ (r3 >> 4)) & 0xf000f000f000f0ULL);
        r3 ^= (r16 | (r16 << 4));
        a[19] = ((r3 >> 32) & 0xffffffffULL);
        r3 = (r3 & 0xffffffffULL);
        a[35] = v_in[7] ^  p_mask[0];
        r16 = ((a[35] ^ (a[35] >> 4))
            & 0xf000f000f000f0ULL);
        a[35] ^= (r16 | (r16 << 4));
        r4 = v_in[22] ^  p_mask[0];
        r16 = ((r4 ^ (r4 >> 4)) & 0xf000f000f000f0ULL);
        r4 ^= (r16 | (r16 << 4));
        a[20] = ((r4 >> 32) & 0xffffffffULL);
        r4 = (r4 & 0xffffffffULL);
        a[36] = v_in[23] ^  p_mask[0];
        r16 = ((a[36] ^ (a[36] >> 4))
            & 0xf000f000f000f0ULL);
        a[36] ^= (r16 | (r16 << 4));
        r5 = v_in[10] ^  p_mask[0];
        r16 = ((r5 ^ (r5 >> 4)) & 0xf000f000f000f0ULL);
        r5 ^= (r16 | (r16 << 4));
        a[21] = ((r5 >> 32) & 0xffffffffULL);
        r5 = (r5 & 0xffffffffULL);
        a[37] = v_in[11] ^  p_mask[0];
        r16 = ((a[37] ^ (a[37] >> 4))
            & 0xf000f000f000f0ULL);
        a[37] ^= (r16 | (r16 << 4));
        r6 = v_in[12] ^  p_mask[0];
        r16 = ((r6 ^ (r6 >> 4)) & 0xf000f000f000f0ULL);
        r6 ^= (r16 | (r16 << 4));
        a[22] = ((r6 >> 32) & 0xffffffffULL);
        r6 = (r6 & 0xffffffffULL);
        a[38] = v_in[13] ^  p_mask[0];
        r16 = ((a[38] ^ (a[38] >> 4))
            & 0xf000f000f000f0ULL);
        a[38] ^= (r16 | (r16 << 4));
        r7 = v_in[16] ^  p_mask[2];
        r16 = ((r7 ^ (r7 >> 4)) & 0xf000f000f000f0ULL);
        r7 ^= (r16 | (r16 << 4));
        a[23] = ((r7 >> 32) & 0xffffffffULL);
        r7 = (r7 & 0xffffffffULL);
        a[39] = v_in[17] ^  p_mask[2];
        r16 = ((a[39] ^ (a[39] >> 4))
            & 0xf000f000f000f0ULL);
        a[39] ^= (r16 | (r16 << 4));
        r8 = v_in[14] ^  p_mask[0];
        r16 = ((r8 ^ (r8 >> 4)) & 0xf000f000f000f0ULL);
        r8 ^= (r16 | (r16 << 4));
        a[24] = ((r8 >> 32) & 0xffffffffULL);
        r8 = (r8 & 0xffffffffULL);
        a[40] = v_in[15] ^  p_mask[0];
        r16 = ((a[40] ^ (a[40] >> 4))
            & 0xf000f000f000f0ULL);
        a[40] ^= (r16 | (r16 << 4));
        r9 = v_in[18] ^  p_mask[0];
        r16 = ((r9 ^ (r9 >> 4)) & 0xf000f000f000f0ULL);
        r9 ^= (r16 | (r16 << 4));
        a[25] = ((r9 >> 32) & 0xffffffffULL);
        r9 = (r9 & 0xffffffffULL);
        a[41] = v_in[19] ^  p_mask[0];
        r16 = ((a[41] ^ (a[41] >> 4))
            & 0xf000f000f000f0ULL);
        a[41] ^= (r16 | (r16 << 4));
        r10 = v_in[20] ^  p_mask[0];
        r16 = ((r10 ^ (r10 >> 4)) & 0xf000f000f000f0ULL);
        r10 ^= (r16 | (r16 << 4));
        a[26] = ((r10 >> 32) & 0xffffffffULL);
        r10 = (r10 & 0xffffffffULL);
        a[42] = v_in[21] ^  p_mask[0];
        r16 = ((a[42] ^ (a[42] >> 4))
            & 0xf000f000f000f0ULL);
        a[42] ^= (r16 | (r16 << 4));
        r11 = v_in[8] ^  p_mask[2];
        r16 = ((r11 ^ (r11 >> 4)) & 0xf000f000f000f0ULL);
        r11 ^= (r16 | (r16 << 4));
        a[27] = ((r11 >> 32) & 0xffffffffULL);
        r11 = (r11 & 0xffffffffULL);
        a[43] = v_in[9] ^  p_mask[2];
        r16 = ((a[43] ^ (a[43] >> 4))
            & 0xf000f000f000f0ULL);
        a[43] ^= (r16 | (r16 << 4));
        r12 = v_in[24] ^  p_mask[0];
        r16 = ((r12 ^ (r12 >> 4)) & 0xf000f000f000f0ULL);
        r12 ^= (r16 | (r16 << 4));
        a[28] = ((r12 >> 32) & 0xffffffffULL);
        r12 = (r12 & 0xffffffffULL);
        a[44] = v_in[25] ^  p_mask[0];
        r16 = ((a[44] ^ (a[44] >> 4))
            & 0xf000f000f000f0ULL);
        a[44] ^= (r16 | (r16 << 4));
        r13 = v_in[4] ^  p_mask[2];
        r16 = ((r13 ^ (r13 >> 4)) & 0xf000f000f000f0ULL);
        r13 ^= (r16 | (r16 << 4));
        a[29] = ((r13 >> 32) & 0xffffffffULL);
        r13 = (r13 & 0xffffffffULL);
        a[45] = v_in[5] ^  p_mask[2];
        r16 = ((a[45] ^ (a[45] >> 4))
            & 0xf000f000f000f0ULL);
        a[45] ^= (r16 | (r16 << 4));
        r14 = v_in[2] ^  p_mask[2];
        r16 = ((r14 ^ (r14 >> 4)) & 0xf000f000f000f0ULL);
        r14 ^= (r16 | (r16 << 4));
        a[30] = ((r14 >> 32) & 0xffffffffULL);
        r14 = (r14 & 0xffffffffULL);
        a[46] = v_in[3] ^  p_mask[2];
        r16 = ((a[46] ^ (a[46] >> 4))
            & 0xf000f000f000f0ULL);
        a[46] ^= (r16 | (r16 << 4));
        r15 = v_in[30] ^  p_mask[2];
        r16 = ((r15 ^ (r15 >> 4)) & 0xf000f000f000f0ULL);
        r15 ^= (r16 | (r16 << 4));
        a[31] = ((r15 >> 32) & 0xffffffffULL);
        r15 = (r15 & 0xffffffffULL);
        a[47] = v_in[31] ^  p_mask[2];
        r16 = ((a[47] ^ (a[47] >> 4))
            & 0xf000f000f000f0ULL);
        a[47] ^= (r16 | (r16 << 4));
        // 128 lines of code, 272 operations
        i = 0;
        goto l_mmv15_op_l64_2;
        l_mmv15_op_l64_1:
        a[(i) + 0] = r0;
        a[(i) + 1] = r1;
        a[(i) + 2] = r2;
        a[(i) + 3] = r3;
        a[(i) + 4] = r4;
        a[(i) + 5] = r5;
        a[(i) + 6] = r6;
        a[(i) + 7] = r7;
        a[(i) + 8] = r8;
        a[(i) + 9] = r9;
        a[(i) + 10] = r10;
        a[(i) + 11] = r11;
        a[(i) + 12] = r12;
        a[(i) + 13] = r13;
        a[(i) + 14] = r14;
        a[(i) + 15] = r15;
        i += 16;
        r0 = a[(i) + 0];
        r1 = a[(i) + 1];
        r2 = a[(i) + 2];
        r3 = a[(i) + 3];
        r4 = a[(i) + 4];
        r5 = a[(i) + 5];
        r6 = a[(i) + 6];
        r7 = a[(i) + 7];
        r8 = a[(i) + 8];
        r9 = a[(i) + 9];
        r10 = a[(i) + 10];
        r11 = a[(i) + 11];
        r12 = a[(i) + 12];
        r13 = a[(i) + 13];
        r14 = a[(i) + 14];
        r15 = a[(i) + 15];
        l_mmv15_op_l64_2:
        // Expansion for Hadamard operation:
        // There is no space for a carry bit between bit fields. So 
        // we move bit field 2*i + 1  to bit field 2*i + 8.
        r0 = ((r0 & 0xf0f0f0fULL)
            | ((r0 & 0xf0f0f0f0ULL) << 28));
        r1 = ((r1 & 0xf0f0f0fULL)
            | ((r1 & 0xf0f0f0f0ULL) << 28));
        r2 = ((r2 & 0xf0f0f0fULL)
            | ((r2 & 0xf0f0f0f0ULL) << 28));
        r3 = ((r3 & 0xf0f0f0fULL)
            | ((r3 & 0xf0f0f0f0ULL) << 28));
        r4 = ((r4 & 0xf0f0f0fULL)
            | ((r4 & 0xf0f0f0f0ULL) << 28));
        r5 = ((r5 & 0xf0f0f0fULL)
            | ((r5 & 0xf0f0f0f0ULL) << 28));
        r6 = ((r6 & 0xf0f0f0fULL)
            | ((r6 & 0xf0f0f0f0ULL) << 28));
        r7 = ((r7 & 0xf0f0f0fULL)
            | ((r7 & 0xf0f0f0f0ULL) << 28));
        r8 = ((r8 & 0xf0f0f0fULL)
            | ((r8 & 0xf0f0f0f0ULL) << 28));
        r9 = ((r9 & 0xf0f0f0fULL)
            | ((r9 & 0xf0f0f0f0ULL) << 28));
        r10 = ((r10 & 0xf0f0f0fULL)
            | ((r10 & 0xf0f0f0f0ULL) << 28));
        r11 = ((r11 & 0xf0f0f0fULL)
            | ((r11 & 0xf0f0f0f0ULL) << 28));
        r12 = ((r12 & 0xf0f0f0fULL)
            | ((r12 & 0xf0f0f0f0ULL) << 28));
        r13 = ((r13 & 0xf0f0f0fULL)
            | ((r13 & 0xf0f0f0f0ULL) << 28));
        r14 = ((r14 & 0xf0f0f0fULL)
            | ((r14 & 0xf0f0f0f0ULL) << 28));
        r15 = ((r15 & 0xf0f0f0fULL)
            | ((r15 & 0xf0f0f0f0ULL) << 28));
        // Vector is now  r(i) for i = 
        //  0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
        // Butterfly: v[i], v[i+2] = v[i]+v[i+2], v[i]-v[i+2]
        r16 = (((r0 << 8) & 0xf000f000f000f00ULL)
            | ((r0 & 0xf000f000f000f00ULL) >> 8));
        r0 = ((r0 ^ 0xf000f000f000f00ULL) + r16);
        r16 = (r0 & 0x1010101010101010ULL);
        r0 = ((r0 - r16) + (r16 >> 4));
        r16 = (((r1 << 8) & 0xf000f000f000f00ULL)
            | ((r1 & 0xf000f000f000f00ULL) >> 8));
        r1 = ((r1 ^ 0xf000f000f000f00ULL) + r16);
        r16 = (r1 & 0x1010101010101010ULL);
        r1 = ((r1 - r16) + (r16 >> 4));
        r16 = (((r2 << 8) & 0xf000f000f000f00ULL)
            | ((r2 & 0xf000f000f000f00ULL) >> 8));
        r2 = ((r2 ^ 0xf000f000f000f00ULL) + r16);
        r16 = (r2 & 0x1010101010101010ULL);
        r2 = ((r2 - r16) + (r16 >> 4));
        r16 = (((r3 << 8) & 0xf000f000f000f00ULL)
            | ((r3 & 0xf000f000f000f00ULL) >> 8));
        r3 = ((r3 ^ 0xf000f000f000f00ULL) + r16);
        r16 = (r3 & 0x1010101010101010ULL);
        r3 = ((r3 - r16) + (r16 >> 4));
        r16 = (((r4 << 8) & 0xf000f000f000f00ULL)
            | ((r4 & 0xf000f000f000f00ULL) >> 8));
        r4 = ((r4 ^ 0xf000f000f000f00ULL) + r16);
        r16 = (r4 & 0x1010101010101010ULL);
        r4 = ((r4 - r16) + (r16 >> 4));
        r16 = (((r5 << 8) & 0xf000f000f000f00ULL)
            | ((r5 & 0xf000f000f000f00ULL) >> 8));
        r5 = ((r5 ^ 0xf000f000f000f00ULL) + r16);
        r16 = (r5 & 0x1010101010101010ULL);
        r5 = ((r5 - r16) + (r16 >> 4));
        r16 = (((r6 << 8) & 0xf000f000f000f00ULL)
            | ((r6 & 0xf000f000f000f00ULL) >> 8));
        r6 = ((r6 ^ 0xf000f000f000f00ULL) + r16);
        r16 = (r6 & 0x1010101010101010ULL);
        r6 = ((r6 - r16) + (r16 >> 4));
        r16 = (((r7 << 8) & 0xf000f000f000f00ULL)
            | ((r7 & 0xf000f000f000f00ULL) >> 8));
        r7 = ((r7 ^ 0xf000f000f000f00ULL) + r16);
        r16 = (r7 & 0x1010101010101010ULL);
        r7 = ((r7 - r16) + (r16 >> 4));
        r16 = (((r8 << 8) & 0xf000f000f000f00ULL)
            | ((r8 & 0xf000f000f000f00ULL) >> 8));
        r8 = ((r8 ^ 0xf000f000f000f00ULL) + r16);
        r16 = (r8 & 0x1010101010101010ULL);
        r8 = ((r8 - r16) + (r16 >> 4));
        r16 = (((r9 << 8) & 0xf000f000f000f00ULL)
            | ((r9 & 0xf000f000f000f00ULL) >> 8));
        r9 = ((r9 ^ 0xf000f000f000f00ULL) + r16);
        r16 = (r9 & 0x1010101010101010ULL);
        r9 = ((r9 - r16) + (r16 >> 4));
        r16 = (((r10 << 8) & 0xf000f000f000f00ULL)
            | ((r10 & 0xf000f000f000f00ULL) >> 8));
        r10 = ((r10 ^ 0xf000f000f000f00ULL) + r16);
        r16 = (r10 & 0x1010101010101010ULL);
        r10 = ((r10 - r16) + (r16 >> 4));
        r16 = (((r11 << 8) & 0xf000f000f000f00ULL)
            | ((r11 & 0xf000f000f000f00ULL) >> 8));
        r11 = ((r11 ^ 0xf000f000f000f00ULL) + r16);
        r16 = (r11 & 0x1010101010101010ULL);
        r11 = ((r11 - r16) + (r16 >> 4));
        r16 = (((r12 << 8) & 0xf000f000f000f00ULL)
            | ((r12 & 0xf000f000f000f00ULL) >> 8));
        r12 = ((r12 ^ 0xf000f000f000f00ULL) + r16);
        r16 = (r12 & 0x1010101010101010ULL);
        r12 = ((r12 - r16) + (r16 >> 4));
        r16 = (((r13 << 8) & 0xf000f000f000f00ULL)
            | ((r13 & 0xf000f000f000f00ULL) >> 8));
        r13 = ((r13 ^ 0xf000f000f000f00ULL) + r16);
        r16 = (r13 & 0x1010101010101010ULL);
        r13 = ((r13 - r16) + (r16 >> 4));
        r16 = (((r14 << 8) & 0xf000f000f000f00ULL)
            | ((r14 & 0xf000f000f000f00ULL) >> 8));
        r14 = ((r14 ^ 0xf000f000f000f00ULL) + r16);
        r16 = (r14 & 0x1010101010101010ULL);
        r14 = ((r14 - r16) + (r16 >> 4));
        r16 = (((r15 << 8) & 0xf000f000f000f00ULL)
            | ((r15 & 0xf000f000f000f00ULL) >> 8));
        r15 = ((r15 ^ 0xf000f000f000f00ULL) + r16);
        r16 = (r15 & 0x1010101010101010ULL);
        r15 = ((r15 - r16) + (r16 >> 4));
        // Vector is now  r(i) for i = 
        //  0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
        // Butterfly: v[i], v[i+8] = v[i]+v[i+8], v[i]-v[i+8]
        r16 = ((r0 << 32) | (r0 >> 32));
        r0 = ((r0 ^ 0xf0f0f0f00000000ULL) + r16);
        r16 = (r0 & 0x1010101010101010ULL);
        r0 = ((r0 - r16) + (r16 >> 4));
        r16 = ((r1 << 32) | (r1 >> 32));
        r1 = ((r1 ^ 0xf0f0f0f00000000ULL) + r16);
        r16 = (r1 & 0x1010101010101010ULL);
        r1 = ((r1 - r16) + (r16 >> 4));
        r16 = ((r2 << 32) | (r2 >> 32));
        r2 = ((r2 ^ 0xf0f0f0f00000000ULL) + r16);
        r16 = (r2 & 0x1010101010101010ULL);
        r2 = ((r2 - r16) + (r16 >> 4));
        r16 = ((r3 << 32) | (r3 >> 32));
        r3 = ((r3 ^ 0xf0f0f0f00000000ULL) + r16);
        r16 = (r3 & 0x1010101010101010ULL);
        r3 = ((r3 - r16) + (r16 >> 4));
        r16 = ((r4 << 32) | (r4 >> 32));
        r4 = ((r4 ^ 0xf0f0f0f00000000ULL) + r16);
        r16 = (r4 & 0x1010101010101010ULL);
        r4 = ((r4 - r16) + (r16 >> 4));
        r16 = ((r5 << 32) | (r5 >> 32));
        r5 = ((r5 ^ 0xf0f0f0f00000000ULL) + r16);
        r16 = (r5 & 0x1010101010101010ULL);
        r5 = ((r5 - r16) + (r16 >> 4));
        r16 = ((r6 << 32) | (r6 >> 32));
        r6 = ((r6 ^ 0xf0f0f0f00000000ULL) + r16);
        r16 = (r6 & 0x1010101010101010ULL);
        r6 = ((r6 - r16) + (r16 >> 4));
        r16 = ((r7 << 32) | (r7 >> 32));
        r7 = ((r7 ^ 0xf0f0f0f00000000ULL) + r16);
        r16 = (r7 & 0x1010101010101010ULL);
        r7 = ((r7 - r16) + (r16 >> 4));
        r16 = ((r8 << 32) | (r8 >> 32));
        r8 = ((r8 ^ 0xf0f0f0f00000000ULL) + r16);
        r16 = (r8 & 0x1010101010101010ULL);
        r8 = ((r8 - r16) + (r16 >> 4));
        r16 = ((r9 << 32) | (r9 >> 32));
        r9 = ((r9 ^ 0xf0f0f0f00000000ULL) + r16);
        r16 = (r9 & 0x1010101010101010ULL);
        r9 = ((r9 - r16) + (r16 >> 4));
        r16 = ((r10 << 32) | (r10 >> 32));
        r10 = ((r10 ^ 0xf0f0f0f00000000ULL) + r16);
        r16 = (r10 & 0x1010101010101010ULL);
        r10 = ((r10 - r16) + (r16 >> 4));
        r16 = ((r11 << 32) | (r11 >> 32));
        r11 = ((r11 ^ 0xf0f0f0f00000000ULL) + r16);
        r16 = (r11 & 0x1010101010101010ULL);
        r11 = ((r11 - r16) + (r16 >> 4));
        r16 = ((r12 << 32) | (r12 >> 32));
        r12 = ((r12 ^ 0xf0f0f0f00000000ULL) + r16);
        r16 = (r12 & 0x1010101010101010ULL);
        r12 = ((r12 - r16) + (r16 >> 4));
        r16 = ((r13 << 32) | (r13 >> 32));
        r13 = ((r13 ^ 0xf0f0f0f00000000ULL) + r16);
        r16 = (r13 & 0x1010101010101010ULL);
        r13 = ((r13 - r16) + (r16 >> 4));
        r16 = ((r14 << 32) | (r14 >> 32));
        r14 = ((r14 ^ 0xf0f0f0f00000000ULL) + r16);
        r16 = (r14 & 0x1010101010101010ULL);
        r14 = ((r14 - r16) + (r16 >> 4));
        r16 = ((r15 << 32) | (r15 >> 32));
        r15 = ((r15 ^ 0xf0f0f0f00000000ULL) + r16);
        r16 = (r15 & 0x1010101010101010ULL);
        r15 = ((r15 - r16) + (r16 >> 4));
        // Vector is now  r(i) for i = 
        //  0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
        // Butterfly: v[i], v[i+16] = v[i]+v[i+16], v[i]-v[i+16]
        r16 = (r0 + (r1 ^ 0xf0f0f0f0f0f0f0fULL));
        r0 = (r0 + r1);
        r1 = (r16 & 0x1010101010101010ULL);
        r1 = ((r16 - r1) + (r1 >> 4));
        r16 = (r0 & 0x1010101010101010ULL);
        r0 = ((r0 - r16) + (r16 >> 4));
        r16 = (r2 + (r3 ^ 0xf0f0f0f0f0f0f0fULL));
        r2 = (r2 + r3);
        r3 = (r16 & 0x1010101010101010ULL);
        r3 = ((r16 - r3) + (r3 >> 4));
        r16 = (r2 & 0x1010101010101010ULL);
        r2 = ((r2 - r16) + (r16 >> 4));
        r16 = (r4 + (r5 ^ 0xf0f0f0f0f0f0f0fULL));
        r4 = (r4 + r5);
        r5 = (r16 & 0x1010101010101010ULL);
        r5 = ((r16 - r5) + (r5 >> 4));
        r16 = (r4 & 0x1010101010101010ULL);
        r4 = ((r4 - r16) + (r16 >> 4));
        r16 = (r6 + (r7 ^ 0xf0f0f0f0f0f0f0fULL));
        r6 = (r6 + r7);
        r7 = (r16 & 0x1010101010101010ULL);
        r7 = ((r16 - r7) + (r7 >> 4));
        r16 = (r6 & 0x1010101010101010ULL);
        r6 = ((r6 - r16) + (r16 >> 4));
        r16 = (r8 + (r9 ^ 0xf0f0f0f0f0f0f0fULL));
        r8 = (r8 + r9);
        r9 = (r16 & 0x1010101010101010ULL);
        r9 = ((r16 - r9) + (r9 >> 4));
        r16 = (r8 & 0x1010101010101010ULL);
        r8 = ((r8 - r16) + (r16 >> 4));
        r16 = (r10 + (r11 ^ 0xf0f0f0f0f0f0f0fULL));
        r10 = (r10 + r11);
        r11 = (r16 & 0x1010101010101010ULL);
        r11 = ((r16 - r11) + (r11 >> 4));
        r16 = (r10 & 0x1010101010101010ULL);
        r10 = ((r10 - r16) + (r16 >> 4));
        r16 = (r12 + (r13 ^ 0xf0f0f0f0f0f0f0fULL));
        r12 = (r12 + r13);
        r13 = (r16 & 0x1010101010101010ULL);
        r13 = ((r16 - r13) + (r13 >> 4));
        r16 = (r12 & 0x1010101010101010ULL);
        r12 = ((r12 - r16) + (r16 >> 4));
        r16 = (r14 + (r15 ^ 0xf0f0f0f0f0f0f0fULL));
        r14 = (r14 + r15);
        r15 = (r16 & 0x1010101010101010ULL);
        r15 = ((r16 - r15) + (r15 >> 4));
        r16 = (r14 & 0x1010101010101010ULL);
        r14 = ((r14 - r16) + (r16 >> 4));
        // Vector is now  r(i) for i = 
        //  0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
        // Butterfly: v[i], v[i+32] = v[i]+v[i+32], v[i]-v[i+32]
        r16 = (r0 + (r2 ^ 0xf0f0f0f0f0f0f0fULL));
        r0 = (r0 + r2);
        r2 = (r16 & 0x1010101010101010ULL);
        r2 = ((r16 - r2) + (r2 >> 4));
        r16 = (r0 & 0x1010101010101010ULL);
        r0 = ((r0 - r16) + (r16 >> 4));
        r16 = (r1 + (r3 ^ 0xf0f0f0f0f0f0f0fULL));
        r1 = (r1 + r3);
        r3 = (r16 & 0x1010101010101010ULL);
        r3 = ((r16 - r3) + (r3 >> 4));
        r16 = (r1 & 0x1010101010101010ULL);
        r1 = ((r1 - r16) + (r16 >> 4));
        r16 = (r4 + (r6 ^ 0xf0f0f0f0f0f0f0fULL));
        r4 = (r4 + r6);
        r6 = (r16 & 0x1010101010101010ULL);
        r6 = ((r16 - r6) + (r6 >> 4));
        r16 = (r4 & 0x1010101010101010ULL);
        r4 = ((r4 - r16) + (r16 >> 4));
        r16 = (r5 + (r7 ^ 0xf0f0f0f0f0f0f0fULL));
        r5 = (r5 + r7);
        r7 = (r16 & 0x1010101010101010ULL);
        r7 = ((r16 - r7) + (r7 >> 4));
        r16 = (r5 & 0x1010101010101010ULL);
        r5 = ((r5 - r16) + (r16 >> 4));
        r16 = (r8 + (r10 ^ 0xf0f0f0f0f0f0f0fULL));
        r8 = (r8 + r10);
        r10 = (r16 & 0x1010101010101010ULL);
        r10 = ((r16 - r10) + (r10 >> 4));
        r16 = (r8 & 0x1010101010101010ULL);
        r8 = ((r8 - r16) + (r16 >> 4));
        r16 = (r9 + (r11 ^ 0xf0f0f0f0f0f0f0fULL));
        r9 = (r9 + r11);
        r11 = (r16 & 0x1010101010101010ULL);
        r11 = ((r16 - r11) + (r11 >> 4));
        r16 = (r9 & 0x1010101010101010ULL);
        r9 = ((r9 - r16) + (r16 >> 4));
        r16 = (r12 + (r14 ^ 0xf0f0f0f0f0f0f0fULL));
        r12 = (r12 + r14);
        r14 = (r16 & 0x1010101010101010ULL);
        r14 = ((r16 - r14) + (r14 >> 4));
        r16 = (r12 & 0x1010101010101010ULL);
        r12 = ((r12 - r16) + (r16 >> 4));
        r16 = (r13 + (r15 ^ 0xf0f0f0f0f0f0f0fULL));
        r13 = (r13 + r15);
        r15 = (r16 & 0x1010101010101010ULL);
        r15 = ((r16 - r15) + (r15 >> 4));
        r16 = (r13 & 0x1010101010101010ULL);
        r13 = ((r13 - r16) + (r16 >> 4));
        // Vector is now  r(i) for i = 
        //  0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
        // Butterfly: v[i], v[i+64] = v[i]+v[i+64], v[i]-v[i+64]
        r16 = (r0 + (r4 ^ 0xf0f0f0f0f0f0f0fULL));
        r0 = (r0 + r4);
        r4 = (r16 & 0x1010101010101010ULL);
        r4 = ((r16 - r4) + (r4 >> 4));
        r16 = (r0 & 0x1010101010101010ULL);
        r0 = ((r0 - r16) + (r16 >> 4));
        r16 = (r1 + (r5 ^ 0xf0f0f0f0f0f0f0fULL));
        r1 = (r1 + r5);
        r5 = (r16 & 0x1010101010101010ULL);
        r5 = ((r16 - r5) + (r5 >> 4));
        r16 = (r1 & 0x1010101010101010ULL);
        r1 = ((r1 - r16) + (r16 >> 4));
        r16 = (r2 + (r6 ^ 0xf0f0f0f0f0f0f0fULL));
        r2 = (r2 + r6);
        r6 = (r16 & 0x1010101010101010ULL);
        r6 = ((r16 - r6) + (r6 >> 4));
        r16 = (r2 & 0x1010101010101010ULL);
        r2 = ((r2 - r16) + (r16 >> 4));
        r16 = (r3 + (r7 ^ 0xf0f0f0f0f0f0f0fULL));
        r3 = (r3 + r7);
        r7 = (r16 & 0x1010101010101010ULL);
        r7 = ((r16 - r7) + (r7 >> 4));
        r16 = (r3 & 0x1010101010101010ULL);
        r3 = ((r3 - r16) + (r16 >> 4));
        r16 = (r8 + (r12 ^ 0xf0f0f0f0f0f0f0fULL));
        r8 = (r8 + r12);
        r12 = (r16 & 0x1010101010101010ULL);
        r12 = ((r16 - r12) + (r12 >> 4));
        r16 = (r8 & 0x1010101010101010ULL);
        r8 = ((r8 - r16) + (r16 >> 4));
        r16 = (r9 + (r13 ^ 0xf0f0f0f0f0f0f0fULL));
        r9 = (r9 + r13);
        r13 = (r16 & 0x1010101010101010ULL);
        r13 = ((r16 - r13) + (r13 >> 4));
        r16 = (r9 & 0x1010101010101010ULL);
        r9 = ((r9 - r16) + (r16 >> 4));
        r16 = (r10 + (r14 ^ 0xf0f0f0f0f0f0f0fULL));
        r10 = (r10 + r14);
        r14 = (r16 & 0x1010101010101010ULL);
        r14 = ((r16 - r14) + (r14 >> 4));
        r16 = (r10 & 0x1010101010101010ULL);
        r10 = ((r10 - r16) + (r16 >> 4));
        r16 = (r11 + (r15 ^ 0xf0f0f0f0f0f0f0fULL));
        r11 = (r11 + r15);
        r15 = (r16 & 0x1010101010101010ULL);
        r15 = ((r16 - r15) + (r15 >> 4));
        r16 = (r11 & 0x1010101010101010ULL);
        r11 = ((r11 - r16) + (r16 >> 4));
        // Vector is now  r(i) for i = 
        //  0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
        // Butterfly: v[i], v[i+128] = v[i]+v[i+128], v[i]-v[i+128]
        r16 = (r0 + (r8 ^ 0xf0f0f0f0f0f0f0fULL));
        r0 = (r0 + r8);
        r8 = (r16 & 0x1010101010101010ULL);
        r8 = ((r16 - r8) + (r8 >> 4));
        r16 = (r0 & 0x1010101010101010ULL);
        r0 = ((r0 - r16) + (r16 >> 4));
        r16 = (r1 + (r9 ^ 0xf0f0f0f0f0f0f0fULL));
        r1 = (r1 + r9);
        r9 = (r16 & 0x1010101010101010ULL);
        r9 = ((r16 - r9) + (r9 >> 4));
        r16 = (r1 & 0x1010101010101010ULL);
        r1 = ((r1 - r16) + (r16 >> 4));
        r16 = (r2 + (r10 ^ 0xf0f0f0f0f0f0f0fULL));
        r2 = (r2 + r10);
        r10 = (r16 & 0x1010101010101010ULL);
        r10 = ((r16 - r10) + (r10 >> 4));
        r16 = (r2 & 0x1010101010101010ULL);
        r2 = ((r2 - r16) + (r16 >> 4));
        r16 = (r3 + (r11 ^ 0xf0f0f0f0f0f0f0fULL));
        r3 = (r3 + r11);
        r11 = (r16 & 0x1010101010101010ULL);
        r11 = ((r16 - r11) + (r11 >> 4));
        r16 = (r3 & 0x1010101010101010ULL);
        r3 = ((r3 - r16) + (r16 >> 4));
        r16 = (r4 + (r12 ^ 0xf0f0f0f0f0f0f0fULL));
        r4 = (r4 + r12);
        r12 = (r16 & 0x1010101010101010ULL);
        r12 = ((r16 - r12) + (r12 >> 4));
        r16 = (r4 & 0x1010101010101010ULL);
        r4 = ((r4 - r16) + (r16 >> 4));
        r16 = (r5 + (r13 ^ 0xf0f0f0f0f0f0f0fULL));
        r5 = (r5 + r13);
        r13 = (r16 & 0x1010101010101010ULL);
        r13 = ((r16 - r13) + (r13 >> 4));
        r16 = (r5 & 0x1010101010101010ULL);
        r5 = ((r5 - r16) + (r16 >> 4));
        r16 = (r6 + (r14 ^ 0xf0f0f0f0f0f0f0fULL));
        r6 = (r6 + r14);
        r14 = (r16 & 0x1010101010101010ULL);
        r14 = ((r16 - r14) + (r14 >> 4));
        r16 = (r6 & 0x1010101010101010ULL);
        r6 = ((r6 - r16) + (r16 >> 4));
        r16 = (r7 + (r15 ^ 0xf0f0f0f0f0f0f0fULL));
        r7 = (r7 + r15);
        r15 = (r16 & 0x1010101010101010ULL);
        r15 = ((r16 - r15) + (r15 >> 4));
        r16 = (r7 & 0x1010101010101010ULL);
        r7 = ((r7 - r16) + (r16 >> 4));
        // Vector is now  r(i) for i = 
        //  0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
        // Reverse expansion for Hadamard operation
        r0 = ((r0 & 0xf0f0f0fULL)
            | ((r0 >> 28) & 0xf0f0f0f0ULL));
        r1 = ((r1 & 0xf0f0f0fULL)
            | ((r1 >> 28) & 0xf0f0f0f0ULL));
        r2 = ((r2 & 0xf0f0f0fULL)
            | ((r2 >> 28) & 0xf0f0f0f0ULL));
        r3 = ((r3 & 0xf0f0f0fULL)
            | ((r3 >> 28) & 0xf0f0f0f0ULL));
        r4 = ((r4 & 0xf0f0f0fULL)
            | ((r4 >> 28) & 0xf0f0f0f0ULL));
        r5 = ((r5 & 0xf0f0f0fULL)
            | ((r5 >> 28) & 0xf0f0f0f0ULL));
        r6 = ((r6 & 0xf0f0f0fULL)
            | ((r6 >> 28) & 0xf0f0f0f0ULL));
        r7 = ((r7 & 0xf0f0f0fULL)
            | ((r7 >> 28) & 0xf0f0f0f0ULL));
        r8 = ((r8 & 0xf0f0f0fULL)
            | ((r8 >> 28) & 0xf0f0f0f0ULL));
        r9 = ((r9 & 0xf0f0f0fULL)
            | ((r9 >> 28) & 0xf0f0f0f0ULL));
        r10 = ((r10 & 0xf0f0f0fULL)
            | ((r10 >> 28) & 0xf0f0f0f0ULL));
        r11 = ((r11 & 0xf0f0f0fULL)
            | ((r11 >> 28) & 0xf0f0f0f0ULL));
        r12 = ((r12 & 0xf0f0f0fULL)
            | ((r12 >> 28) & 0xf0f0f0f0ULL));
        r13 = ((r13 & 0xf0f0f0fULL)
            | ((r13 >> 28) & 0xf0f0f0f0ULL));
        r14 = ((r14 & 0xf0f0f0fULL)
            | ((r14 >> 28) & 0xf0f0f0f0ULL));
        r15 = ((r15 & 0xf0f0f0fULL)
            | ((r15 >> 28) & 0xf0f0f0f0ULL));
        // Vector is now  r(i) for i = 
        //  0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
        // 352 lines of code, 800 operations
        if (i < 32) goto l_mmv15_op_l64_1;
        a[0] = ((a[0] & 0xffffffffULL)
            | ((a[16] << 32) & 0xffffffff00000000ULL));
        r0 = r0;
        a[0] = (((a[0] & 0x7777777777777777ULL) << 1)
            | ((a[0] & 0x8888888888888888ULL) >> 3));
        r0 = (((r0 & 0x7777777777777777ULL) << 1)
            | ((r0 & 0x8888888888888888ULL) >> 3));
        v_out[0] = a[0] ^  p_mask[6];
        v_out[1] = (r0 ^  p_mask[6]) & 0xffffffffULL;
        a[1] = ((a[1] & 0xffffffffULL)
            | ((a[17] << 32) & 0xffffffff00000000ULL));
        r1 = r1;
        a[1] = (((a[1] & 0x7777777777777777ULL) << 1)
            | ((a[1] & 0x8888888888888888ULL) >> 3));
        r1 = (((r1 & 0x7777777777777777ULL) << 1)
            | ((r1 & 0x8888888888888888ULL) >> 3));
        v_out[2] = a[1] ^  p_mask[6];
        v_out[3] = (r1 ^  p_mask[6]) & 0xffffffffULL;
        a[2] = ((a[2] & 0xffffffffULL)
            | ((a[18] << 32) & 0xffffffff00000000ULL));
        r2 = r2;
        a[2] = (((a[2] & 0x7777777777777777ULL) << 1)
            | ((a[2] & 0x8888888888888888ULL) >> 3));
        r2 = (((r2 & 0x7777777777777777ULL) << 1)
            | ((r2 & 0x8888888888888888ULL) >> 3));
        v_out[4] = a[2] ^  p_mask[6];
        v_out[5] = (r2 ^  p_mask[6]) & 0xffffffffULL;
        a[3] = ((a[3] & 0xffffffffULL)
            | ((a[19] << 32) & 0xffffffff00000000ULL));
        r3 = r3;
        a[3] = (((a[3] & 0x7777777777777777ULL) << 1)
            | ((a[3] & 0x8888888888888888ULL) >> 3));
        r3 = (((r3 & 0x7777777777777777ULL) << 1)
            | ((r3 & 0x8888888888888888ULL) >> 3));
        v_out[6] = a[3] ^  p_mask[4];
        v_out[7] = (r3 ^  p_mask[4]) & 0xffffffffULL;
        a[4] = ((a[4] & 0xffffffffULL)
            | ((a[20] << 32) & 0xffffffff00000000ULL));
        r4 = r4;
        a[4] = (((a[4] & 0x7777777777777777ULL) << 1)
            | ((a[4] & 0x8888888888888888ULL) >> 3));
        r4 = (((r4 & 0x7777777777777777ULL) << 1)
            | ((r4 & 0x8888888888888888ULL) >> 3));
        v_out[8] = a[4] ^  p_mask[6];
        v_out[9] = (r4 ^  p_mask[6]) & 0xffffffffULL;
        a[5] = ((a[5] & 0xffffffffULL)
            | ((a[21] << 32) & 0xffffffff00000000ULL));
        r5 = r5;
        a[5] = (((a[5] & 0x7777777777777777ULL) << 1)
            | ((a[5] & 0x8888888888888888ULL) >> 3));
        r5 = (((r5 & 0x7777777777777777ULL) << 1)
            | ((r5 & 0x8888888888888888ULL) >> 3));
        v_out[10] = a[5] ^  p_mask[4];
        v_out[11] = (r5 ^  p_mask[4]) & 0xffffffffULL;
        a[6] = ((a[6] & 0xffffffffULL)
            | ((a[22] << 32) & 0xffffffff00000000ULL));
        r6 = r6;
        a[6] = (((a[6] & 0x7777777777777777ULL) << 1)
            | ((a[6] & 0x8888888888888888ULL) >> 3));
        r6 = (((r6 & 0x7777777777777777ULL) << 1)
            | ((r6 & 0x8888888888888888ULL) >> 3));
        v_out[12] = a[6] ^  p_mask[4];
        v_out[13] = (r6 ^  p_mask[4]) & 0xffffffffULL;
        a[7] = ((a[7] & 0xffffffffULL)
            | ((a[23] << 32) & 0xffffffff00000000ULL));
        r7 = r7;
        a[7] = (((a[7] & 0x7777777777777777ULL) << 1)
            | ((a[7] & 0x8888888888888888ULL) >> 3));
        r7 = (((r7 & 0x7777777777777777ULL) << 1)
            | ((r7 & 0x8888888888888888ULL) >> 3));
        v_out[14] = a[7] ^  p_mask[4];
        v_out[15] = (r7 ^  p_mask[4]) & 0xffffffffULL;
        a[8] = ((a[8] & 0xffffffffULL)
            | ((a[24] << 32) & 0xffffffff00000000ULL));
        r8 = r8;
        a[8] = (((a[8] & 0x7777777777777777ULL) << 1)
            | ((a[8] & 0x8888888888888888ULL) >> 3));
        r8 = (((r8 & 0x7777777777777777ULL) << 1)
            | ((r8 & 0x8888888888888888ULL) >> 3));
        v_out[16] = a[8] ^  p_mask[6];
        v_out[17] = (r8 ^  p_mask[6]) & 0xffffffffULL;
        a[9] = ((a[9] & 0xffffffffULL)
            | ((a[25] << 32) & 0xffffffff00000000ULL));
        r9 = r9;
        a[9] = (((a[9] & 0x7777777777777777ULL) << 1)
            | ((a[9] & 0x8888888888888888ULL) >> 3));
        r9 = (((r9 & 0x7777777777777777ULL) << 1)
            | ((r9 & 0x8888888888888888ULL) >> 3));
        v_out[18] = a[9] ^  p_mask[4];
        v_out[19] = (r9 ^  p_mask[4]) & 0xffffffffULL;
        a[10] = ((a[10] & 0xffffffffULL)
            | ((a[26] << 32) & 0xffffffff00000000ULL));
        r10 = r10;
        a[10] = (((a[10] & 0x7777777777777777ULL) << 1)
            | ((a[10] & 0x8888888888888888ULL) >> 3));
        r10 = (((r10 & 0x7777777777777777ULL) << 1)
            | ((r10 & 0x8888888888888888ULL) >> 3));
        v_out[20] = a[10] ^  p_mask[4];
        v_out[21] = (r10 ^  p_mask[4]) & 0xffffffffULL;
        a[11] = ((a[11] & 0xffffffffULL)
            | ((a[27] << 32) & 0xffffffff00000000ULL));
        r11 = r11;
        a[11] = (((a[11] & 0x7777777777777777ULL) << 1)
            | ((a[11] & 0x8888888888888888ULL) >> 3));
        r11 = (((r11 & 0x7777777777777777ULL) << 1)
            | ((r11 & 0x8888888888888888ULL) >> 3));
        v_out[22] = a[11] ^  p_mask[4];
        v_out[23] = (r11 ^  p_mask[4]) & 0xffffffffULL;
        a[12] = ((a[12] & 0xffffffffULL)
            | ((a[28] << 32) & 0xffffffff00000000ULL));
        r12 = r12;
        a[12] = (((a[12] & 0x7777777777777777ULL) << 1)
            | ((a[12] & 0x8888888888888888ULL) >> 3));
        r12 = (((r12 & 0x7777777777777777ULL) << 1)
            | ((r12 & 0x8888888888888888ULL) >> 3));
        v_out[24] = a[12] ^  p_mask[4];
        v_out[25] = (r12 ^  p_mask[4]) & 0xffffffffULL;
        a[13] = ((a[13] & 0xffffffffULL)
            | ((a[29] << 32) & 0xffffffff00000000ULL));
        r13 = r13;
        a[13] = (((a[13] & 0x7777777777777777ULL) << 1)
            | ((a[13] & 0x8888888888888888ULL) >> 3));
        r13 = (((r13 & 0x7777777777777777ULL) << 1)
            | ((r13 & 0x8888888888888888ULL) >> 3));
        v_out[26] = a[13] ^  p_mask[4];
        v_out[27] = (r13 ^  p_mask[4]) & 0xffffffffULL;
        a[14] = ((a[14] & 0xffffffffULL)
            | ((a[30] << 32) & 0xffffffff00000000ULL));
        r14 = r14;
        a[14] = (((a[14] & 0x7777777777777777ULL) << 1)
            | ((a[14] & 0x8888888888888888ULL) >> 3));
        r14 = (((r14 & 0x7777777777777777ULL) << 1)
            | ((r14 & 0x8888888888888888ULL) >> 3));
        v_out[28] = a[14] ^  p_mask[4];
        v_out[29] = (r14 ^  p_mask[4]) & 0xffffffffULL;
        a[15] = ((a[15] & 0xffffffffULL)
            | ((a[31] << 32) & 0xffffffff00000000ULL));
        r15 = r15;
        a[15] = (((a[15] & 0x7777777777777777ULL) << 1)
            | ((a[15] & 0x8888888888888888ULL) >> 3));
        r15 = (((r15 & 0x7777777777777777ULL) << 1)
            | ((r15 & 0x8888888888888888ULL) >> 3));
        v_out[30] = a[15] ^  p_mask[6];
        v_out[31] = (r15 ^  p_mask[6]) & 0xffffffffULL;
        v_in += 32;
        v_out += 32;
        // 96 lines of code, 272 operations
        }
        // End of automatically generated matrix operation.
 
    }
}  


static void mm_op15_xi_a(uint_mmv_t *v_in,  uint32_t exp1, uint_mmv_t *v_out)
{
    uint_fast32_t i;
    uint_mmv_t e_mask =  -((uint_mmv_t)exp1 & 0x1ULL);
    for (i = 0; i < 6; ++i) {
        // %%MUL_MATRIX_XI16 v_in, e_mask, v_out

        // This is an automatically generated matrix operation, do not change!
        {
        uint_mmv_t r0, r1, r2, r3, r4;
        uint_mmv_t r5, r6, r7, r8;

        uint_fast32_t i;
        // TODO: write comment!!!
        // 
        for (i = 0; i < 2; ++i) {
        e_mask = ~(e_mask);
        r0 = v_in[0] ^  (0xfff0fff0fff0fff0ULL & e_mask);
        r4 = ((r0 ^ (r0 >> 4)) & 0xf000f000f000f0ULL);
        r0 ^= (r4 | (r4 << 4));
        r1 = v_in[4] ^  (0xf000f000f000fULL & e_mask);
        r4 = ((r1 ^ (r1 >> 4)) & 0xf000f000f000f0ULL);
        r1 ^= (r4 | (r4 << 4));
        r2 = v_in[2] ^  (0xf000f000f000fULL & e_mask);
        r4 = ((r2 ^ (r2 >> 4)) & 0xf000f000f000f0ULL);
        r2 ^= (r4 | (r4 << 4));
        r3 = v_in[6] ^  (0xf000f000f000fULL & e_mask);
        r4 = ((r3 ^ (r3 >> 4)) & 0xf000f000f000f0ULL);
        r3 ^= (r4 | (r4 << 4));
        // Expansion for Hadamard operation:
        // There is no space for a carry bit between bit fields. So 
        // we move bit field 2*i + 1  to bit field 2*i + 64.
        r4 = ((r0 >> 4) & 0xf0f0f0f0f0f0f0fULL);
        r0 = (r0 & 0xf0f0f0f0f0f0f0fULL);
        r5 = ((r1 >> 4) & 0xf0f0f0f0f0f0f0fULL);
        r1 = (r1 & 0xf0f0f0f0f0f0f0fULL);
        r6 = ((r2 >> 4) & 0xf0f0f0f0f0f0f0fULL);
        r2 = (r2 & 0xf0f0f0f0f0f0f0fULL);
        r7 = ((r3 >> 4) & 0xf0f0f0f0f0f0f0fULL);
        r3 = (r3 & 0xf0f0f0f0f0f0f0fULL);
        // Vector is now  r(i) for i = 0,1,2,3,4,5,6,7
        // Butterfly: v[i], v[i+2] = v[i]+v[i+2], v[i]-v[i+2]
        r8 = (((r0 << 8) & 0xf000f000f000f00ULL)
            | ((r0 & 0xf000f000f000f00ULL) >> 8));
        r0 = ((r0 ^ 0xf000f000f000f00ULL) + r8);
        r8 = (r0 & 0x1010101010101010ULL);
        r0 = ((r0 - r8) + (r8 >> 4));
        r8 = (((r1 << 8) & 0xf000f000f000f00ULL)
            | ((r1 & 0xf000f000f000f00ULL) >> 8));
        r1 = ((r1 ^ 0xf000f000f000f00ULL) + r8);
        r8 = (r1 & 0x1010101010101010ULL);
        r1 = ((r1 - r8) + (r8 >> 4));
        r8 = (((r2 << 8) & 0xf000f000f000f00ULL)
            | ((r2 & 0xf000f000f000f00ULL) >> 8));
        r2 = ((r2 ^ 0xf000f000f000f00ULL) + r8);
        r8 = (r2 & 0x1010101010101010ULL);
        r2 = ((r2 - r8) + (r8 >> 4));
        r8 = (((r3 << 8) & 0xf000f000f000f00ULL)
            | ((r3 & 0xf000f000f000f00ULL) >> 8));
        r3 = ((r3 ^ 0xf000f000f000f00ULL) + r8);
        r8 = (r3 & 0x1010101010101010ULL);
        r3 = ((r3 - r8) + (r8 >> 4));
        r8 = (((r4 << 8) & 0xf000f000f000f00ULL)
            | ((r4 & 0xf000f000f000f00ULL) >> 8));
        r4 = ((r4 ^ 0xf000f000f000f00ULL) + r8);
        r8 = (r4 & 0x1010101010101010ULL);
        r4 = ((r4 - r8) + (r8 >> 4));
        r8 = (((r5 << 8) & 0xf000f000f000f00ULL)
            | ((r5 & 0xf000f000f000f00ULL) >> 8));
        r5 = ((r5 ^ 0xf000f000f000f00ULL) + r8);
        r8 = (r5 & 0x1010101010101010ULL);
        r5 = ((r5 - r8) + (r8 >> 4));
        r8 = (((r6 << 8) & 0xf000f000f000f00ULL)
            | ((r6 & 0xf000f000f000f00ULL) >> 8));
        r6 = ((r6 ^ 0xf000f000f000f00ULL) + r8);
        r8 = (r6 & 0x1010101010101010ULL);
        r6 = ((r6 - r8) + (r8 >> 4));
        r8 = (((r7 << 8) & 0xf000f000f000f00ULL)
            | ((r7 & 0xf000f000f000f00ULL) >> 8));
        r7 = ((r7 ^ 0xf000f000f000f00ULL) + r8);
        r8 = (r7 & 0x1010101010101010ULL);
        r7 = ((r7 - r8) + (r8 >> 4));
        // Vector is now  r(i) for i = 0,1,2,3,4,5,6,7
        // Butterfly: v[i], v[i+16] = v[i]+v[i+16], v[i]-v[i+16]
        r8 = (r0 + (r1 ^ 0xf0f0f0f0f0f0f0fULL));
        r0 = (r0 + r1);
        r1 = (r8 & 0x1010101010101010ULL);
        r1 = ((r8 - r1) + (r1 >> 4));
        r8 = (r0 & 0x1010101010101010ULL);
        r0 = ((r0 - r8) + (r8 >> 4));
        r8 = (r2 + (r3 ^ 0xf0f0f0f0f0f0f0fULL));
        r2 = (r2 + r3);
        r3 = (r8 & 0x1010101010101010ULL);
        r3 = ((r8 - r3) + (r3 >> 4));
        r8 = (r2 & 0x1010101010101010ULL);
        r2 = ((r2 - r8) + (r8 >> 4));
        r8 = (r4 + (r5 ^ 0xf0f0f0f0f0f0f0fULL));
        r4 = (r4 + r5);
        r5 = (r8 & 0x1010101010101010ULL);
        r5 = ((r8 - r5) + (r5 >> 4));
        r8 = (r4 & 0x1010101010101010ULL);
        r4 = ((r4 - r8) + (r8 >> 4));
        r8 = (r6 + (r7 ^ 0xf0f0f0f0f0f0f0fULL));
        r6 = (r6 + r7);
        r7 = (r8 & 0x1010101010101010ULL);
        r7 = ((r8 - r7) + (r7 >> 4));
        r8 = (r6 & 0x1010101010101010ULL);
        r6 = ((r6 - r8) + (r8 >> 4));
        // Vector is now  r(i) for i = 0,1,2,3,4,5,6,7
        // Butterfly: v[i], v[i+32] = v[i]+v[i+32], v[i]-v[i+32]
        r8 = (r0 + (r2 ^ 0xf0f0f0f0f0f0f0fULL));
        r0 = (r0 + r2);
        r2 = (r8 & 0x1010101010101010ULL);
        r2 = ((r8 - r2) + (r2 >> 4));
        r8 = (r0 & 0x1010101010101010ULL);
        r0 = ((r0 - r8) + (r8 >> 4));
        r8 = (r1 + (r3 ^ 0xf0f0f0f0f0f0f0fULL));
        r1 = (r1 + r3);
        r3 = (r8 & 0x1010101010101010ULL);
        r3 = ((r8 - r3) + (r3 >> 4));
        r8 = (r1 & 0x1010101010101010ULL);
        r1 = ((r1 - r8) + (r8 >> 4));
        r8 = (r4 + (r6 ^ 0xf0f0f0f0f0f0f0fULL));
        r4 = (r4 + r6);
        r6 = (r8 & 0x1010101010101010ULL);
        r6 = ((r8 - r6) + (r6 >> 4));
        r8 = (r4 & 0x1010101010101010ULL);
        r4 = ((r4 - r8) + (r8 >> 4));
        r8 = (r5 + (r7 ^ 0xf0f0f0f0f0f0f0fULL));
        r5 = (r5 + r7);
        r7 = (r8 & 0x1010101010101010ULL);
        r7 = ((r8 - r7) + (r7 >> 4));
        r8 = (r5 & 0x1010101010101010ULL);
        r5 = ((r5 - r8) + (r8 >> 4));
        // Vector is now  r(i) for i = 0,1,2,3,4,5,6,7
        // Butterfly: v[i], v[i+64] = v[i]+v[i+64], v[i]-v[i+64]
        r8 = (r0 + (r4 ^ 0xf0f0f0f0f0f0f0fULL));
        r0 = (r0 + r4);
        r4 = (r8 & 0x1010101010101010ULL);
        r4 = ((r8 - r4) + (r4 >> 4));
        r8 = (r0 & 0x1010101010101010ULL);
        r0 = ((r0 - r8) + (r8 >> 4));
        r8 = (r1 + (r5 ^ 0xf0f0f0f0f0f0f0fULL));
        r1 = (r1 + r5);
        r5 = (r8 & 0x1010101010101010ULL);
        r5 = ((r8 - r5) + (r5 >> 4));
        r8 = (r1 & 0x1010101010101010ULL);
        r1 = ((r1 - r8) + (r8 >> 4));
        r8 = (r2 + (r6 ^ 0xf0f0f0f0f0f0f0fULL));
        r2 = (r2 + r6);
        r6 = (r8 & 0x1010101010101010ULL);
        r6 = ((r8 - r6) + (r6 >> 4));
        r8 = (r2 & 0x1010101010101010ULL);
        r2 = ((r2 - r8) + (r8 >> 4));
        r8 = (r3 + (r7 ^ 0xf0f0f0f0f0f0f0fULL));
        r3 = (r3 + r7);
        r7 = (r8 & 0x1010101010101010ULL);
        r7 = ((r8 - r7) + (r7 >> 4));
        r8 = (r3 & 0x1010101010101010ULL);
        r3 = ((r3 - r8) + (r8 >> 4));
        // Vector is now  r(i) for i = 0,1,2,3,4,5,6,7
        // Reverse expansion for Hadamard operation
        r0 ^= (r4 << 4);
        r1 ^= (r5 << 4);
        r2 ^= (r6 << 4);
        r3 ^= (r7 << 4);
        // Vector is now  r(i) for i = 0,1,2,3
        e_mask = ~(e_mask);
        r0 = (((r0 & 0x3333333333333333ULL) << 2)
            | ((r0 & 0xccccccccccccccccULL) >> 2));
        v_out[0] = r0 ^ (e_mask & 0xfff0fff0fff0fff0ULL);
        r1 = (((r1 & 0x3333333333333333ULL) << 2)
            | ((r1 & 0xccccccccccccccccULL) >> 2));
        v_out[2] = r1 ^ (e_mask & 0xf000f000f000fULL);
        r2 = (((r2 & 0x3333333333333333ULL) << 2)
            | ((r2 & 0xccccccccccccccccULL) >> 2));
        v_out[4] = r2 ^ (e_mask & 0xf000f000f000fULL);
        r3 = (((r3 & 0x3333333333333333ULL) << 2)
            | ((r3 & 0xccccccccccccccccULL) >> 2));
        v_out[6] = r3 ^ (e_mask & 0xf000f000f000fULL);
        // 138 lines of code, 302 operations
        v_in++;
        v_out++;
        }
        v_out[-1] &= 0xffffffffULL;
        v_out[1] &= 0xffffffffULL;
        v_out[3] &= 0xffffffffULL;
        v_out[5] &= 0xffffffffULL;
        v_in += 6;
        v_out += 6;
        }
        // End of automatically generated matrix operation.
 
    }
}  


// %%EXPORT p
void mm_op15_xi(uint_mmv_t *v_in,  uint32_t exp, uint_mmv_t *v_out)
{
    uint_mmv_t i, exp1;
 
    exp %= 3;
    if (exp == 0) {
        for (i = 0; i < 15468; ++i) v_out[i] = v_in[i];
        return;
    }
    exp1 =  (uint_mmv_t)exp - 0x1ULL;

    // Do monomial part, i.e. tags B, C, T, X
    mm_op15_xi_mon(v_in, exp1, v_out);

    // Do tag A
    mm_op15_xi_a(v_in, exp1, v_out); 

    // Do tags X, Y
    for (i = 0; i < 4; ++i) {
        uint_mmv_t *p_src = v_in + MM_OP15_OFS_Z + (i << HALF_YZ_SHIFT);
        mm_op15_xi_yz(p_src, exp1, v_out + TAB15_XI64_OFFSET[exp1][i]);
    }
    
}


