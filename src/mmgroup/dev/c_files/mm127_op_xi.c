/////////////////////////////////////////////////////////////////////////////
// This C file has been created automatically. Do not edit!!!
/////////////////////////////////////////////////////////////////////////////

// %%COMMENT
// TODO: Yet to be documented!!!


#include <string.h>
#include "mat24_functions.h"
#include "mm_op127.h"   


static void mm_op127_xi_mon(uint_mmv_t *v_in,  uint32_t exp1, uint_mmv_t *v_out)
{
    uint_mmv_t *p_src, *p_dest;
    uint_fast32_t i, j;
    uint_fast32_t diff = exp1 ? 4096 : 0;
    uint8_t b[2496], *p_b;
    mm_sub_table_xi_type *p_tables = mm_sub_table_xi[exp1];
    uint16_t *p_perm;
    uint32_t *p_sign;



    ///////////////////////////////////////////////////////////////
    // Map tag BC to tag BC.
    ///////////////////////////////////////////////////////////////
    p_src = v_in + 96;
    p_dest = v_out + 96;
    p_sign = p_tables[0].p_sign;
    p_perm = p_tables[0].p_perm;

    for (i = 0; i < 1; ++i) {
        p_b = b;
        for (j = 0; j < 78; ++j) {
           // %%OP_XI_LOAD p_src, p_b, {i_src.SHAPE[2]}
           uint_mmv_t r0;
           r0 =  (p_src)[0];
           r0 &= 0x7f7f7f7f7f7f7f7fULL;
           (p_b)[0] = r0 >> 0;
           (p_b)[1] = r0 >> 8;
           (p_b)[2] = r0 >> 16;
           (p_b)[3] = r0 >> 24;
           (p_b)[4] = r0 >> 32;
           (p_b)[5] = r0 >> 40;
           (p_b)[6] = r0 >> 48;
           (p_b)[7] = r0 >> 56;
           r0 =  (p_src)[1];
           r0 &= 0x7f7f7f7f7f7f7f7fULL;
           (p_b)[8] = r0 >> 0;
           (p_b)[9] = r0 >> 8;
           (p_b)[10] = r0 >> 16;
           (p_b)[11] = r0 >> 24;
           (p_b)[12] = r0 >> 32;
           (p_b)[13] = r0 >> 40;
           (p_b)[14] = r0 >> 48;
           (p_b)[15] = r0 >> 56;
           r0 =  (p_src)[2];
           r0 &= 0x7f7f7f7f7f7f7f7fULL;
           (p_b)[16] = r0 >> 0;
           (p_b)[17] = r0 >> 8;
           (p_b)[18] = r0 >> 16;
           (p_b)[19] = r0 >> 24;
           (p_b)[20] = r0 >> 32;
           (p_b)[21] = r0 >> 40;
           (p_b)[22] = r0 >> 48;
           (p_b)[23] = r0 >> 56;
           r0 =  (p_src)[3];
           r0 &= 0x7f7f7f7f7f7f7f7fULL;
           (p_b)[24] = r0 >> 0;
           (p_b)[25] = r0 >> 8;
           (p_b)[26] = r0 >> 16;
           (p_b)[27] = r0 >> 24;
           (p_b)[28] = r0 >> 32;
           (p_b)[29] = r0 >> 40;
           (p_b)[30] = r0 >> 48;
           (p_b)[31] = r0 >> 56;
           p_src += 4;
           p_b += 32;
        }

        for (j = 0; j < 78; ++j) {
           // %%OP_XI_STORE b, p_perm, p_sign, p_dest, {i_dest.SHAPE[2]}
           uint_mmv_t r0, r1;
           r0 = ((uint_mmv_t)(b[p_perm[0]]) << 0)
             + ((uint_mmv_t)(b[p_perm[1]]) << 8)
             + ((uint_mmv_t)(b[p_perm[2]]) << 16)
             + ((uint_mmv_t)(b[p_perm[3]]) << 24)
             + ((uint_mmv_t)(b[p_perm[4]]) << 32)
             + ((uint_mmv_t)(b[p_perm[5]]) << 40)
             + ((uint_mmv_t)(b[p_perm[6]]) << 48)
             + ((uint_mmv_t)(b[p_perm[7]]) << 56);
           r1 = (p_sign)[0] >> 0;
           // Spread bits 0,...,7 of r1 to the (8-bit long) fields
           // of r1. A field of r1 is set to 0x7f if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xfULL) + ((r1 & 0xf0ULL) << 28);
           r1 = (r1 & 0x300000003ULL) 
               +  ((r1 & 0xc0000000cULL) << 14);
           r1 = (r1 & 0x1000100010001ULL) 
               +  ((r1 & 0x2000200020002ULL) << 7);
           r1 = (((r1) << 7) - (r1));
           // Bit spreading done.
           (p_dest)[0] = r0 ^ r1;
           r0 = ((uint_mmv_t)(b[p_perm[8]]) << 0)
             + ((uint_mmv_t)(b[p_perm[9]]) << 8)
             + ((uint_mmv_t)(b[p_perm[10]]) << 16)
             + ((uint_mmv_t)(b[p_perm[11]]) << 24)
             + ((uint_mmv_t)(b[p_perm[12]]) << 32)
             + ((uint_mmv_t)(b[p_perm[13]]) << 40)
             + ((uint_mmv_t)(b[p_perm[14]]) << 48)
             + ((uint_mmv_t)(b[p_perm[15]]) << 56);
           r1 = (p_sign)[0] >> 8;
           // Spread bits 0,...,7 of r1 to the (8-bit long) fields
           // of r1. A field of r1 is set to 0x7f if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xfULL) + ((r1 & 0xf0ULL) << 28);
           r1 = (r1 & 0x300000003ULL) 
               +  ((r1 & 0xc0000000cULL) << 14);
           r1 = (r1 & 0x1000100010001ULL) 
               +  ((r1 & 0x2000200020002ULL) << 7);
           r1 = (((r1) << 7) - (r1));
           // Bit spreading done.
           (p_dest)[1] = r0 ^ r1;
           r0 = ((uint_mmv_t)(b[p_perm[16]]) << 0)
             + ((uint_mmv_t)(b[p_perm[17]]) << 8)
             + ((uint_mmv_t)(b[p_perm[18]]) << 16)
             + ((uint_mmv_t)(b[p_perm[19]]) << 24)
             + ((uint_mmv_t)(b[p_perm[20]]) << 32)
             + ((uint_mmv_t)(b[p_perm[21]]) << 40)
             + ((uint_mmv_t)(b[p_perm[22]]) << 48)
             + ((uint_mmv_t)(b[p_perm[23]]) << 56);
           r1 = (p_sign)[0] >> 16;
           // Spread bits 0,...,7 of r1 to the (8-bit long) fields
           // of r1. A field of r1 is set to 0x7f if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xfULL) + ((r1 & 0xf0ULL) << 28);
           r1 = (r1 & 0x300000003ULL) 
               +  ((r1 & 0xc0000000cULL) << 14);
           r1 = (r1 & 0x1000100010001ULL) 
               +  ((r1 & 0x2000200020002ULL) << 7);
           r1 = (((r1) << 7) - (r1));
           // Bit spreading done.
           (p_dest)[2] = r0 ^ r1;
           r0 = ((uint_mmv_t)(b[p_perm[24]]) << 0)
             + ((uint_mmv_t)(b[p_perm[25]]) << 8)
             + ((uint_mmv_t)(b[p_perm[26]]) << 16)
             + ((uint_mmv_t)(b[p_perm[27]]) << 24)
             + ((uint_mmv_t)(b[p_perm[28]]) << 32)
             + ((uint_mmv_t)(b[p_perm[29]]) << 40)
             + ((uint_mmv_t)(b[p_perm[30]]) << 48)
             + ((uint_mmv_t)(b[p_perm[31]]) << 56);
           r1 = (p_sign)[0] >> 24;
           // Spread bits 0,...,7 of r1 to the (8-bit long) fields
           // of r1. A field of r1 is set to 0x7f if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xfULL) + ((r1 & 0xf0ULL) << 28);
           r1 = (r1 & 0x300000003ULL) 
               +  ((r1 & 0xc0000000cULL) << 14);
           r1 = (r1 & 0x1000100010001ULL) 
               +  ((r1 & 0x2000200020002ULL) << 7);
           r1 = (((r1) << 7) - (r1));
           // Bit spreading done.
           (p_dest)[3] = r0 ^ r1;
           p_dest += 4;
           p_perm += 32;
           p_sign += 1;
        }
    }

    ///////////////////////////////////////////////////////////////
    // Map tag T0 to tag T0.
    ///////////////////////////////////////////////////////////////
    p_src = v_in + 408;
    p_dest = v_out + 408;
    p_sign = p_tables[1].p_sign;
    p_perm = p_tables[1].p_perm;

    for (i = 0; i < 45; ++i) {
        p_b = b;
        for (j = 0; j < 16; ++j) {
           // %%OP_XI_LOAD p_src, p_b, {i_src.SHAPE[2]}
           uint_mmv_t r0;
           r0 =  (p_src)[0];
           r0 &= 0x7f7f7f7f7f7f7f7fULL;
           (p_b)[0] = r0 >> 0;
           (p_b)[1] = r0 >> 8;
           (p_b)[2] = r0 >> 16;
           (p_b)[3] = r0 >> 24;
           (p_b)[4] = r0 >> 32;
           (p_b)[5] = r0 >> 40;
           (p_b)[6] = r0 >> 48;
           (p_b)[7] = r0 >> 56;
           r0 =  (p_src)[1];
           r0 &= 0x7f7f7f7f7f7f7f7fULL;
           (p_b)[8] = r0 >> 0;
           (p_b)[9] = r0 >> 8;
           (p_b)[10] = r0 >> 16;
           (p_b)[11] = r0 >> 24;
           (p_b)[12] = r0 >> 32;
           (p_b)[13] = r0 >> 40;
           (p_b)[14] = r0 >> 48;
           (p_b)[15] = r0 >> 56;
           r0 =  (p_src)[2];
           r0 &= 0x7f7f7f7f7f7f7f7fULL;
           (p_b)[16] = r0 >> 0;
           (p_b)[17] = r0 >> 8;
           (p_b)[18] = r0 >> 16;
           (p_b)[19] = r0 >> 24;
           (p_b)[20] = r0 >> 32;
           (p_b)[21] = r0 >> 40;
           (p_b)[22] = r0 >> 48;
           (p_b)[23] = r0 >> 56;
           r0 =  (p_src)[3];
           r0 &= 0x7f7f7f7f7f7f7f7fULL;
           (p_b)[24] = r0 >> 0;
           (p_b)[25] = r0 >> 8;
           (p_b)[26] = r0 >> 16;
           (p_b)[27] = r0 >> 24;
           (p_b)[28] = r0 >> 32;
           (p_b)[29] = r0 >> 40;
           (p_b)[30] = r0 >> 48;
           (p_b)[31] = r0 >> 56;
           p_src += 4;
           p_b += 32;
        }

        for (j = 0; j < 16; ++j) {
           // %%OP_XI_STORE b, p_perm, p_sign, p_dest, {i_dest.SHAPE[2]}
           uint_mmv_t r0, r1;
           r0 = ((uint_mmv_t)(b[p_perm[0]]) << 0)
             + ((uint_mmv_t)(b[p_perm[1]]) << 8)
             + ((uint_mmv_t)(b[p_perm[2]]) << 16)
             + ((uint_mmv_t)(b[p_perm[3]]) << 24)
             + ((uint_mmv_t)(b[p_perm[4]]) << 32)
             + ((uint_mmv_t)(b[p_perm[5]]) << 40)
             + ((uint_mmv_t)(b[p_perm[6]]) << 48)
             + ((uint_mmv_t)(b[p_perm[7]]) << 56);
           r1 = (p_sign)[0] >> 0;
           // Spread bits 0,...,7 of r1 to the (8-bit long) fields
           // of r1. A field of r1 is set to 0x7f if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xfULL) + ((r1 & 0xf0ULL) << 28);
           r1 = (r1 & 0x300000003ULL) 
               +  ((r1 & 0xc0000000cULL) << 14);
           r1 = (r1 & 0x1000100010001ULL) 
               +  ((r1 & 0x2000200020002ULL) << 7);
           r1 = (((r1) << 7) - (r1));
           // Bit spreading done.
           (p_dest)[0] = r0 ^ r1;
           r0 = ((uint_mmv_t)(b[p_perm[8]]) << 0)
             + ((uint_mmv_t)(b[p_perm[9]]) << 8)
             + ((uint_mmv_t)(b[p_perm[10]]) << 16)
             + ((uint_mmv_t)(b[p_perm[11]]) << 24)
             + ((uint_mmv_t)(b[p_perm[12]]) << 32)
             + ((uint_mmv_t)(b[p_perm[13]]) << 40)
             + ((uint_mmv_t)(b[p_perm[14]]) << 48)
             + ((uint_mmv_t)(b[p_perm[15]]) << 56);
           r1 = (p_sign)[0] >> 8;
           // Spread bits 0,...,7 of r1 to the (8-bit long) fields
           // of r1. A field of r1 is set to 0x7f if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xfULL) + ((r1 & 0xf0ULL) << 28);
           r1 = (r1 & 0x300000003ULL) 
               +  ((r1 & 0xc0000000cULL) << 14);
           r1 = (r1 & 0x1000100010001ULL) 
               +  ((r1 & 0x2000200020002ULL) << 7);
           r1 = (((r1) << 7) - (r1));
           // Bit spreading done.
           (p_dest)[1] = r0 ^ r1;
           r0 = ((uint_mmv_t)(b[p_perm[16]]) << 0)
             + ((uint_mmv_t)(b[p_perm[17]]) << 8)
             + ((uint_mmv_t)(b[p_perm[18]]) << 16)
             + ((uint_mmv_t)(b[p_perm[19]]) << 24)
             + ((uint_mmv_t)(b[p_perm[20]]) << 32)
             + ((uint_mmv_t)(b[p_perm[21]]) << 40)
             + ((uint_mmv_t)(b[p_perm[22]]) << 48)
             + ((uint_mmv_t)(b[p_perm[23]]) << 56);
           r1 = (p_sign)[0] >> 16;
           // Spread bits 0,...,7 of r1 to the (8-bit long) fields
           // of r1. A field of r1 is set to 0x7f if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xfULL) + ((r1 & 0xf0ULL) << 28);
           r1 = (r1 & 0x300000003ULL) 
               +  ((r1 & 0xc0000000cULL) << 14);
           r1 = (r1 & 0x1000100010001ULL) 
               +  ((r1 & 0x2000200020002ULL) << 7);
           r1 = (((r1) << 7) - (r1));
           // Bit spreading done.
           (p_dest)[2] = r0 ^ r1;
           r0 = ((uint_mmv_t)(b[p_perm[24]]) << 0)
             + ((uint_mmv_t)(b[p_perm[25]]) << 8)
             + ((uint_mmv_t)(b[p_perm[26]]) << 16)
             + ((uint_mmv_t)(b[p_perm[27]]) << 24)
             + ((uint_mmv_t)(b[p_perm[28]]) << 32)
             + ((uint_mmv_t)(b[p_perm[29]]) << 40)
             + ((uint_mmv_t)(b[p_perm[30]]) << 48)
             + ((uint_mmv_t)(b[p_perm[31]]) << 56);
           r1 = (p_sign)[0] >> 24;
           // Spread bits 0,...,7 of r1 to the (8-bit long) fields
           // of r1. A field of r1 is set to 0x7f if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xfULL) + ((r1 & 0xf0ULL) << 28);
           r1 = (r1 & 0x300000003ULL) 
               +  ((r1 & 0xc0000000cULL) << 14);
           r1 = (r1 & 0x1000100010001ULL) 
               +  ((r1 & 0x2000200020002ULL) << 7);
           r1 = (((r1) << 7) - (r1));
           // Bit spreading done.
           (p_dest)[3] = r0 ^ r1;
           p_dest += 4;
           p_perm += 32;
           p_sign += 1;
        }
    }

    ///////////////////////////////////////////////////////////////
    // Map tag T1 to tag X0 if e = 1
    // Map tag T1 to tag X1 if e = 2
    ///////////////////////////////////////////////////////////////
    p_src = v_in + 3288;
    p_dest = v_out + 6360;
    p_dest += diff;
    p_sign = p_tables[2].p_sign;
    p_perm = p_tables[2].p_perm;

    for (i = 0; i < 64; ++i) {
        p_b = b;
        for (j = 0; j < 12; ++j) {
           // %%OP_XI_LOAD p_src, p_b, {i_src.SHAPE[2]}
           uint_mmv_t r0;
           r0 =  (p_src)[0];
           r0 &= 0x7f7f7f7f7f7f7f7fULL;
           (p_b)[0] = r0 >> 0;
           (p_b)[1] = r0 >> 8;
           (p_b)[2] = r0 >> 16;
           (p_b)[3] = r0 >> 24;
           (p_b)[4] = r0 >> 32;
           (p_b)[5] = r0 >> 40;
           (p_b)[6] = r0 >> 48;
           (p_b)[7] = r0 >> 56;
           r0 =  (p_src)[1];
           r0 &= 0x7f7f7f7f7f7f7f7fULL;
           (p_b)[8] = r0 >> 0;
           (p_b)[9] = r0 >> 8;
           (p_b)[10] = r0 >> 16;
           (p_b)[11] = r0 >> 24;
           (p_b)[12] = r0 >> 32;
           (p_b)[13] = r0 >> 40;
           (p_b)[14] = r0 >> 48;
           (p_b)[15] = r0 >> 56;
           r0 =  (p_src)[2];
           r0 &= 0x7f7f7f7f7f7f7f7fULL;
           (p_b)[16] = r0 >> 0;
           (p_b)[17] = r0 >> 8;
           (p_b)[18] = r0 >> 16;
           (p_b)[19] = r0 >> 24;
           (p_b)[20] = r0 >> 32;
           (p_b)[21] = r0 >> 40;
           (p_b)[22] = r0 >> 48;
           (p_b)[23] = r0 >> 56;
           r0 =  (p_src)[3];
           r0 &= 0x7f7f7f7f7f7f7f7fULL;
           (p_b)[24] = r0 >> 0;
           (p_b)[25] = r0 >> 8;
           (p_b)[26] = r0 >> 16;
           (p_b)[27] = r0 >> 24;
           (p_b)[28] = r0 >> 32;
           (p_b)[29] = r0 >> 40;
           (p_b)[30] = r0 >> 48;
           (p_b)[31] = r0 >> 56;
           p_src += 4;
           p_b += 32;
        }

        for (j = 0; j < 16; ++j) {
           // %%OP_XI_STORE b, p_perm, p_sign, p_dest, {i_dest.SHAPE[2]}
           uint_mmv_t r0, r1;
           r0 = ((uint_mmv_t)(b[p_perm[0]]) << 0)
             + ((uint_mmv_t)(b[p_perm[1]]) << 8)
             + ((uint_mmv_t)(b[p_perm[2]]) << 16)
             + ((uint_mmv_t)(b[p_perm[3]]) << 24)
             + ((uint_mmv_t)(b[p_perm[4]]) << 32)
             + ((uint_mmv_t)(b[p_perm[5]]) << 40)
             + ((uint_mmv_t)(b[p_perm[6]]) << 48)
             + ((uint_mmv_t)(b[p_perm[7]]) << 56);
           r1 = (p_sign)[0] >> 0;
           // Spread bits 0,...,7 of r1 to the (8-bit long) fields
           // of r1. A field of r1 is set to 0x7f if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xfULL) + ((r1 & 0xf0ULL) << 28);
           r1 = (r1 & 0x300000003ULL) 
               +  ((r1 & 0xc0000000cULL) << 14);
           r1 = (r1 & 0x1000100010001ULL) 
               +  ((r1 & 0x2000200020002ULL) << 7);
           r1 = (((r1) << 7) - (r1));
           // Bit spreading done.
           (p_dest)[0] = r0 ^ r1;
           r0 = ((uint_mmv_t)(b[p_perm[8]]) << 0)
             + ((uint_mmv_t)(b[p_perm[9]]) << 8)
             + ((uint_mmv_t)(b[p_perm[10]]) << 16)
             + ((uint_mmv_t)(b[p_perm[11]]) << 24)
             + ((uint_mmv_t)(b[p_perm[12]]) << 32)
             + ((uint_mmv_t)(b[p_perm[13]]) << 40)
             + ((uint_mmv_t)(b[p_perm[14]]) << 48)
             + ((uint_mmv_t)(b[p_perm[15]]) << 56);
           r1 = (p_sign)[0] >> 8;
           // Spread bits 0,...,7 of r1 to the (8-bit long) fields
           // of r1. A field of r1 is set to 0x7f if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xfULL) + ((r1 & 0xf0ULL) << 28);
           r1 = (r1 & 0x300000003ULL) 
               +  ((r1 & 0xc0000000cULL) << 14);
           r1 = (r1 & 0x1000100010001ULL) 
               +  ((r1 & 0x2000200020002ULL) << 7);
           r1 = (((r1) << 7) - (r1));
           // Bit spreading done.
           (p_dest)[1] = r0 ^ r1;
           r0 = ((uint_mmv_t)(b[p_perm[16]]) << 0)
             + ((uint_mmv_t)(b[p_perm[17]]) << 8)
             + ((uint_mmv_t)(b[p_perm[18]]) << 16)
             + ((uint_mmv_t)(b[p_perm[19]]) << 24)
             + ((uint_mmv_t)(b[p_perm[20]]) << 32)
             + ((uint_mmv_t)(b[p_perm[21]]) << 40)
             + ((uint_mmv_t)(b[p_perm[22]]) << 48)
             + ((uint_mmv_t)(b[p_perm[23]]) << 56);
           r1 = (p_sign)[0] >> 16;
           // Spread bits 0,...,7 of r1 to the (8-bit long) fields
           // of r1. A field of r1 is set to 0x7f if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xfULL) + ((r1 & 0xf0ULL) << 28);
           r1 = (r1 & 0x300000003ULL) 
               +  ((r1 & 0xc0000000cULL) << 14);
           r1 = (r1 & 0x1000100010001ULL) 
               +  ((r1 & 0x2000200020002ULL) << 7);
           r1 = (((r1) << 7) - (r1));
           // Bit spreading done.
           (p_dest)[2] = r0 ^ r1;
           p_dest[3] = 0;
           p_dest += 4;
           p_perm += 24;
           p_sign += 1;
        }
    }

    ///////////////////////////////////////////////////////////////
    // Map tag X0 to tag X1 if e = 1
    // Map tag X1 to tag X0 if e = 2
    ///////////////////////////////////////////////////////////////
    p_src = v_in + 6360;
    p_src += diff;
    p_dest = v_out + 10456;
    p_dest -= diff;
    p_sign = p_tables[3].p_sign;
    p_perm = p_tables[3].p_perm;

    for (i = 0; i < 64; ++i) {
        p_b = b;
        for (j = 0; j < 16; ++j) {
           // %%OP_XI_LOAD p_src, p_b, {i_src.SHAPE[2]}
           uint_mmv_t r0;
           r0 =  (p_src)[0];
           r0 &= 0x7f7f7f7f7f7f7f7fULL;
           (p_b)[0] = r0 >> 0;
           (p_b)[1] = r0 >> 8;
           (p_b)[2] = r0 >> 16;
           (p_b)[3] = r0 >> 24;
           (p_b)[4] = r0 >> 32;
           (p_b)[5] = r0 >> 40;
           (p_b)[6] = r0 >> 48;
           (p_b)[7] = r0 >> 56;
           r0 =  (p_src)[1];
           r0 &= 0x7f7f7f7f7f7f7f7fULL;
           (p_b)[8] = r0 >> 0;
           (p_b)[9] = r0 >> 8;
           (p_b)[10] = r0 >> 16;
           (p_b)[11] = r0 >> 24;
           (p_b)[12] = r0 >> 32;
           (p_b)[13] = r0 >> 40;
           (p_b)[14] = r0 >> 48;
           (p_b)[15] = r0 >> 56;
           r0 =  (p_src)[2];
           r0 &= 0x7f7f7f7f7f7f7f7fULL;
           (p_b)[16] = r0 >> 0;
           (p_b)[17] = r0 >> 8;
           (p_b)[18] = r0 >> 16;
           (p_b)[19] = r0 >> 24;
           (p_b)[20] = r0 >> 32;
           (p_b)[21] = r0 >> 40;
           (p_b)[22] = r0 >> 48;
           (p_b)[23] = r0 >> 56;
           p_src += 4;
           p_b += 32;
        }

        for (j = 0; j < 16; ++j) {
           // %%OP_XI_STORE b, p_perm, p_sign, p_dest, {i_dest.SHAPE[2]}
           uint_mmv_t r0, r1;
           r0 = ((uint_mmv_t)(b[p_perm[0]]) << 0)
             + ((uint_mmv_t)(b[p_perm[1]]) << 8)
             + ((uint_mmv_t)(b[p_perm[2]]) << 16)
             + ((uint_mmv_t)(b[p_perm[3]]) << 24)
             + ((uint_mmv_t)(b[p_perm[4]]) << 32)
             + ((uint_mmv_t)(b[p_perm[5]]) << 40)
             + ((uint_mmv_t)(b[p_perm[6]]) << 48)
             + ((uint_mmv_t)(b[p_perm[7]]) << 56);
           r1 = (p_sign)[0] >> 0;
           // Spread bits 0,...,7 of r1 to the (8-bit long) fields
           // of r1. A field of r1 is set to 0x7f if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xfULL) + ((r1 & 0xf0ULL) << 28);
           r1 = (r1 & 0x300000003ULL) 
               +  ((r1 & 0xc0000000cULL) << 14);
           r1 = (r1 & 0x1000100010001ULL) 
               +  ((r1 & 0x2000200020002ULL) << 7);
           r1 = (((r1) << 7) - (r1));
           // Bit spreading done.
           (p_dest)[0] = r0 ^ r1;
           r0 = ((uint_mmv_t)(b[p_perm[8]]) << 0)
             + ((uint_mmv_t)(b[p_perm[9]]) << 8)
             + ((uint_mmv_t)(b[p_perm[10]]) << 16)
             + ((uint_mmv_t)(b[p_perm[11]]) << 24)
             + ((uint_mmv_t)(b[p_perm[12]]) << 32)
             + ((uint_mmv_t)(b[p_perm[13]]) << 40)
             + ((uint_mmv_t)(b[p_perm[14]]) << 48)
             + ((uint_mmv_t)(b[p_perm[15]]) << 56);
           r1 = (p_sign)[0] >> 8;
           // Spread bits 0,...,7 of r1 to the (8-bit long) fields
           // of r1. A field of r1 is set to 0x7f if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xfULL) + ((r1 & 0xf0ULL) << 28);
           r1 = (r1 & 0x300000003ULL) 
               +  ((r1 & 0xc0000000cULL) << 14);
           r1 = (r1 & 0x1000100010001ULL) 
               +  ((r1 & 0x2000200020002ULL) << 7);
           r1 = (((r1) << 7) - (r1));
           // Bit spreading done.
           (p_dest)[1] = r0 ^ r1;
           r0 = ((uint_mmv_t)(b[p_perm[16]]) << 0)
             + ((uint_mmv_t)(b[p_perm[17]]) << 8)
             + ((uint_mmv_t)(b[p_perm[18]]) << 16)
             + ((uint_mmv_t)(b[p_perm[19]]) << 24)
             + ((uint_mmv_t)(b[p_perm[20]]) << 32)
             + ((uint_mmv_t)(b[p_perm[21]]) << 40)
             + ((uint_mmv_t)(b[p_perm[22]]) << 48)
             + ((uint_mmv_t)(b[p_perm[23]]) << 56);
           r1 = (p_sign)[0] >> 16;
           // Spread bits 0,...,7 of r1 to the (8-bit long) fields
           // of r1. A field of r1 is set to 0x7f if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xfULL) + ((r1 & 0xf0ULL) << 28);
           r1 = (r1 & 0x300000003ULL) 
               +  ((r1 & 0xc0000000cULL) << 14);
           r1 = (r1 & 0x1000100010001ULL) 
               +  ((r1 & 0x2000200020002ULL) << 7);
           r1 = (((r1) << 7) - (r1));
           // Bit spreading done.
           (p_dest)[2] = r0 ^ r1;
           p_dest[3] = 0;
           p_dest += 4;
           p_perm += 24;
           p_sign += 1;
        }
    }

    ///////////////////////////////////////////////////////////////
    // Map tag X1 to tag T1 if e = 1
    // Map tag X0 to tag T1 if e = 2
    ///////////////////////////////////////////////////////////////
    p_src = v_in + 10456;
    p_src -= diff;
    p_dest = v_out + 3288;
    p_sign = p_tables[4].p_sign;
    p_perm = p_tables[4].p_perm;

    for (i = 0; i < 64; ++i) {
        p_b = b;
        for (j = 0; j < 16; ++j) {
           // %%OP_XI_LOAD p_src, p_b, {i_src.SHAPE[2]}
           uint_mmv_t r0;
           r0 =  (p_src)[0];
           r0 &= 0x7f7f7f7f7f7f7f7fULL;
           (p_b)[0] = r0 >> 0;
           (p_b)[1] = r0 >> 8;
           (p_b)[2] = r0 >> 16;
           (p_b)[3] = r0 >> 24;
           (p_b)[4] = r0 >> 32;
           (p_b)[5] = r0 >> 40;
           (p_b)[6] = r0 >> 48;
           (p_b)[7] = r0 >> 56;
           r0 =  (p_src)[1];
           r0 &= 0x7f7f7f7f7f7f7f7fULL;
           (p_b)[8] = r0 >> 0;
           (p_b)[9] = r0 >> 8;
           (p_b)[10] = r0 >> 16;
           (p_b)[11] = r0 >> 24;
           (p_b)[12] = r0 >> 32;
           (p_b)[13] = r0 >> 40;
           (p_b)[14] = r0 >> 48;
           (p_b)[15] = r0 >> 56;
           r0 =  (p_src)[2];
           r0 &= 0x7f7f7f7f7f7f7f7fULL;
           (p_b)[16] = r0 >> 0;
           (p_b)[17] = r0 >> 8;
           (p_b)[18] = r0 >> 16;
           (p_b)[19] = r0 >> 24;
           (p_b)[20] = r0 >> 32;
           (p_b)[21] = r0 >> 40;
           (p_b)[22] = r0 >> 48;
           (p_b)[23] = r0 >> 56;
           p_src += 4;
           p_b += 32;
        }

        for (j = 0; j < 12; ++j) {
           // %%OP_XI_STORE b, p_perm, p_sign, p_dest, {i_dest.SHAPE[2]}
           uint_mmv_t r0, r1;
           r0 = ((uint_mmv_t)(b[p_perm[0]]) << 0)
             + ((uint_mmv_t)(b[p_perm[1]]) << 8)
             + ((uint_mmv_t)(b[p_perm[2]]) << 16)
             + ((uint_mmv_t)(b[p_perm[3]]) << 24)
             + ((uint_mmv_t)(b[p_perm[4]]) << 32)
             + ((uint_mmv_t)(b[p_perm[5]]) << 40)
             + ((uint_mmv_t)(b[p_perm[6]]) << 48)
             + ((uint_mmv_t)(b[p_perm[7]]) << 56);
           r1 = (p_sign)[0] >> 0;
           // Spread bits 0,...,7 of r1 to the (8-bit long) fields
           // of r1. A field of r1 is set to 0x7f if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xfULL) + ((r1 & 0xf0ULL) << 28);
           r1 = (r1 & 0x300000003ULL) 
               +  ((r1 & 0xc0000000cULL) << 14);
           r1 = (r1 & 0x1000100010001ULL) 
               +  ((r1 & 0x2000200020002ULL) << 7);
           r1 = (((r1) << 7) - (r1));
           // Bit spreading done.
           (p_dest)[0] = r0 ^ r1;
           r0 = ((uint_mmv_t)(b[p_perm[8]]) << 0)
             + ((uint_mmv_t)(b[p_perm[9]]) << 8)
             + ((uint_mmv_t)(b[p_perm[10]]) << 16)
             + ((uint_mmv_t)(b[p_perm[11]]) << 24)
             + ((uint_mmv_t)(b[p_perm[12]]) << 32)
             + ((uint_mmv_t)(b[p_perm[13]]) << 40)
             + ((uint_mmv_t)(b[p_perm[14]]) << 48)
             + ((uint_mmv_t)(b[p_perm[15]]) << 56);
           r1 = (p_sign)[0] >> 8;
           // Spread bits 0,...,7 of r1 to the (8-bit long) fields
           // of r1. A field of r1 is set to 0x7f if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xfULL) + ((r1 & 0xf0ULL) << 28);
           r1 = (r1 & 0x300000003ULL) 
               +  ((r1 & 0xc0000000cULL) << 14);
           r1 = (r1 & 0x1000100010001ULL) 
               +  ((r1 & 0x2000200020002ULL) << 7);
           r1 = (((r1) << 7) - (r1));
           // Bit spreading done.
           (p_dest)[1] = r0 ^ r1;
           r0 = ((uint_mmv_t)(b[p_perm[16]]) << 0)
             + ((uint_mmv_t)(b[p_perm[17]]) << 8)
             + ((uint_mmv_t)(b[p_perm[18]]) << 16)
             + ((uint_mmv_t)(b[p_perm[19]]) << 24)
             + ((uint_mmv_t)(b[p_perm[20]]) << 32)
             + ((uint_mmv_t)(b[p_perm[21]]) << 40)
             + ((uint_mmv_t)(b[p_perm[22]]) << 48)
             + ((uint_mmv_t)(b[p_perm[23]]) << 56);
           r1 = (p_sign)[0] >> 16;
           // Spread bits 0,...,7 of r1 to the (8-bit long) fields
           // of r1. A field of r1 is set to 0x7f if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xfULL) + ((r1 & 0xf0ULL) << 28);
           r1 = (r1 & 0x300000003ULL) 
               +  ((r1 & 0xc0000000cULL) << 14);
           r1 = (r1 & 0x1000100010001ULL) 
               +  ((r1 & 0x2000200020002ULL) << 7);
           r1 = (((r1) << 7) - (r1));
           // Bit spreading done.
           (p_dest)[2] = r0 ^ r1;
           r0 = ((uint_mmv_t)(b[p_perm[24]]) << 0)
             + ((uint_mmv_t)(b[p_perm[25]]) << 8)
             + ((uint_mmv_t)(b[p_perm[26]]) << 16)
             + ((uint_mmv_t)(b[p_perm[27]]) << 24)
             + ((uint_mmv_t)(b[p_perm[28]]) << 32)
             + ((uint_mmv_t)(b[p_perm[29]]) << 40)
             + ((uint_mmv_t)(b[p_perm[30]]) << 48)
             + ((uint_mmv_t)(b[p_perm[31]]) << 56);
           r1 = (p_sign)[0] >> 24;
           // Spread bits 0,...,7 of r1 to the (8-bit long) fields
           // of r1. A field of r1 is set to 0x7f if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xfULL) + ((r1 & 0xf0ULL) << 28);
           r1 = (r1 & 0x300000003ULL) 
               +  ((r1 & 0xc0000000cULL) << 14);
           r1 = (r1 & 0x1000100010001ULL) 
               +  ((r1 & 0x2000200020002ULL) << 7);
           r1 = (((r1) << 7) - (r1));
           // Bit spreading done.
           (p_dest)[3] = r0 ^ r1;
           p_dest += 4;
           p_perm += 32;
           p_sign += 1;
        }
    }
}

static uint_mmv_t TAB127_XI64_MASK[] = {
// %%TABLE TABLE_MUL_MATRIX_XI64, uint{INT_BITS}
0x0000007f0000007fULL,0x0000000000000000ULL,
0x0000007f0000007fULL,0x7f7f7f7f7f7f7f7fULL,
0x0000000000000000ULL,0x0000007f0000007fULL,
0x7f7f7f7f7f7f7f7fULL,0x0000007f0000007fULL
};


#define HALF_YZ_SHIFT 12

static uint32_t TAB127_XI64_OFFSET[2][4] = {
    {
        MM_OP127_OFS_Z + (2 << HALF_YZ_SHIFT),
        MM_OP127_OFS_Z + (0 << HALF_YZ_SHIFT),
        MM_OP127_OFS_Z + (1 << HALF_YZ_SHIFT),
        MM_OP127_OFS_Z + (3 << HALF_YZ_SHIFT)
    },
    {
        MM_OP127_OFS_Z + (1 << HALF_YZ_SHIFT),
        MM_OP127_OFS_Z + (2 << HALF_YZ_SHIFT),
        MM_OP127_OFS_Z + (0 << HALF_YZ_SHIFT),
        MM_OP127_OFS_Z + (3 << HALF_YZ_SHIFT)
    },
};




static void mm_op127_xi_yz(uint_mmv_t *v_in,  uint32_t exp1, uint_mmv_t *v_out)
{
    uint_fast32_t i;
    uint_mmv_t *p_mask =  TAB127_XI64_MASK + exp1;
    for (i = 0; i < 64; ++i) {
        // %%MUL_MATRIX_XI64 v_in, p_mask, v_out

        // This is an automatically generated matrix operation, do not change!
        {
        uint_mmv_t r0, r1, r2, r3, r4;
        uint_mmv_t r5, r6, r7, r8, r9;
        uint_mmv_t r10, r11, r12, r13, r14;
        uint_mmv_t r15, r16;

        uint_fast32_t i;

        // TODO: write comment!!!
        // 
        for (i = 0; i < 3; ++i) {
        r0 = v_in[0] ^  p_mask[2];
        r16 = ((r0 ^ (r0 >> 8)) & 0x7f0000007f00ULL);
        r0 ^= (r16 | (r16 << 8));
        r1 = v_in[56] ^  p_mask[0];
        r16 = ((r1 ^ (r1 >> 8)) & 0x7f0000007f00ULL);
        r1 ^= (r16 | (r16 << 8));
        r2 = v_in[52] ^  p_mask[0];
        r16 = ((r2 ^ (r2 >> 8)) & 0x7f0000007f00ULL);
        r2 ^= (r16 | (r16 << 8));
        r3 = v_in[12] ^  p_mask[0];
        r16 = ((r3 ^ (r3 >> 8)) & 0x7f0000007f00ULL);
        r3 ^= (r16 | (r16 << 8));
        r4 = v_in[44] ^  p_mask[0];
        r16 = ((r4 ^ (r4 >> 8)) & 0x7f0000007f00ULL);
        r4 ^= (r16 | (r16 << 8));
        r5 = v_in[20] ^  p_mask[0];
        r16 = ((r5 ^ (r5 >> 8)) & 0x7f0000007f00ULL);
        r5 ^= (r16 | (r16 << 8));
        r6 = v_in[24] ^  p_mask[0];
        r16 = ((r6 ^ (r6 >> 8)) & 0x7f0000007f00ULL);
        r6 ^= (r16 | (r16 << 8));
        r7 = v_in[32] ^  p_mask[2];
        r16 = ((r7 ^ (r7 >> 8)) & 0x7f0000007f00ULL);
        r7 ^= (r16 | (r16 << 8));
        r8 = v_in[28] ^  p_mask[0];
        r16 = ((r8 ^ (r8 >> 8)) & 0x7f0000007f00ULL);
        r8 ^= (r16 | (r16 << 8));
        r9 = v_in[36] ^  p_mask[0];
        r16 = ((r9 ^ (r9 >> 8)) & 0x7f0000007f00ULL);
        r9 ^= (r16 | (r16 << 8));
        r10 = v_in[40] ^  p_mask[0];
        r16 = ((r10 ^ (r10 >> 8)) & 0x7f0000007f00ULL);
        r10 ^= (r16 | (r16 << 8));
        r11 = v_in[16] ^  p_mask[2];
        r16 = ((r11 ^ (r11 >> 8)) & 0x7f0000007f00ULL);
        r11 ^= (r16 | (r16 << 8));
        r12 = v_in[48] ^  p_mask[0];
        r16 = ((r12 ^ (r12 >> 8)) & 0x7f0000007f00ULL);
        r12 ^= (r16 | (r16 << 8));
        r13 = v_in[8] ^  p_mask[2];
        r16 = ((r13 ^ (r13 >> 8)) & 0x7f0000007f00ULL);
        r13 ^= (r16 | (r16 << 8));
        r14 = v_in[4] ^  p_mask[2];
        r16 = ((r14 ^ (r14 >> 8)) & 0x7f0000007f00ULL);
        r14 ^= (r16 | (r16 << 8));
        r15 = v_in[60] ^  p_mask[2];
        r16 = ((r15 ^ (r15 >> 8)) & 0x7f0000007f00ULL);
        r15 ^= (r16 | (r16 << 8));
        // Butterfly: v[i], v[i+1] = v[i]+v[i+1], v[i]-v[i+1]
        r16 = (((r0 << 8) & 0x7f007f007f007f00ULL)
            | ((r0 & 0x7f007f007f007f00ULL) >> 8));
        r0 = ((r0 ^ 0x7f007f007f007f00ULL) + r16);
        r16 = (r0 & 0x8080808080808080ULL);
        r0 = ((r0 - r16) + (r16 >> 7));
        r16 = (((r1 << 8) & 0x7f007f007f007f00ULL)
            | ((r1 & 0x7f007f007f007f00ULL) >> 8));
        r1 = ((r1 ^ 0x7f007f007f007f00ULL) + r16);
        r16 = (r1 & 0x8080808080808080ULL);
        r1 = ((r1 - r16) + (r16 >> 7));
        r16 = (((r2 << 8) & 0x7f007f007f007f00ULL)
            | ((r2 & 0x7f007f007f007f00ULL) >> 8));
        r2 = ((r2 ^ 0x7f007f007f007f00ULL) + r16);
        r16 = (r2 & 0x8080808080808080ULL);
        r2 = ((r2 - r16) + (r16 >> 7));
        r16 = (((r3 << 8) & 0x7f007f007f007f00ULL)
            | ((r3 & 0x7f007f007f007f00ULL) >> 8));
        r3 = ((r3 ^ 0x7f007f007f007f00ULL) + r16);
        r16 = (r3 & 0x8080808080808080ULL);
        r3 = ((r3 - r16) + (r16 >> 7));
        r16 = (((r4 << 8) & 0x7f007f007f007f00ULL)
            | ((r4 & 0x7f007f007f007f00ULL) >> 8));
        r4 = ((r4 ^ 0x7f007f007f007f00ULL) + r16);
        r16 = (r4 & 0x8080808080808080ULL);
        r4 = ((r4 - r16) + (r16 >> 7));
        r16 = (((r5 << 8) & 0x7f007f007f007f00ULL)
            | ((r5 & 0x7f007f007f007f00ULL) >> 8));
        r5 = ((r5 ^ 0x7f007f007f007f00ULL) + r16);
        r16 = (r5 & 0x8080808080808080ULL);
        r5 = ((r5 - r16) + (r16 >> 7));
        r16 = (((r6 << 8) & 0x7f007f007f007f00ULL)
            | ((r6 & 0x7f007f007f007f00ULL) >> 8));
        r6 = ((r6 ^ 0x7f007f007f007f00ULL) + r16);
        r16 = (r6 & 0x8080808080808080ULL);
        r6 = ((r6 - r16) + (r16 >> 7));
        r16 = (((r7 << 8) & 0x7f007f007f007f00ULL)
            | ((r7 & 0x7f007f007f007f00ULL) >> 8));
        r7 = ((r7 ^ 0x7f007f007f007f00ULL) + r16);
        r16 = (r7 & 0x8080808080808080ULL);
        r7 = ((r7 - r16) + (r16 >> 7));
        r16 = (((r8 << 8) & 0x7f007f007f007f00ULL)
            | ((r8 & 0x7f007f007f007f00ULL) >> 8));
        r8 = ((r8 ^ 0x7f007f007f007f00ULL) + r16);
        r16 = (r8 & 0x8080808080808080ULL);
        r8 = ((r8 - r16) + (r16 >> 7));
        r16 = (((r9 << 8) & 0x7f007f007f007f00ULL)
            | ((r9 & 0x7f007f007f007f00ULL) >> 8));
        r9 = ((r9 ^ 0x7f007f007f007f00ULL) + r16);
        r16 = (r9 & 0x8080808080808080ULL);
        r9 = ((r9 - r16) + (r16 >> 7));
        r16 = (((r10 << 8) & 0x7f007f007f007f00ULL)
            | ((r10 & 0x7f007f007f007f00ULL) >> 8));
        r10 = ((r10 ^ 0x7f007f007f007f00ULL) + r16);
        r16 = (r10 & 0x8080808080808080ULL);
        r10 = ((r10 - r16) + (r16 >> 7));
        r16 = (((r11 << 8) & 0x7f007f007f007f00ULL)
            | ((r11 & 0x7f007f007f007f00ULL) >> 8));
        r11 = ((r11 ^ 0x7f007f007f007f00ULL) + r16);
        r16 = (r11 & 0x8080808080808080ULL);
        r11 = ((r11 - r16) + (r16 >> 7));
        r16 = (((r12 << 8) & 0x7f007f007f007f00ULL)
            | ((r12 & 0x7f007f007f007f00ULL) >> 8));
        r12 = ((r12 ^ 0x7f007f007f007f00ULL) + r16);
        r16 = (r12 & 0x8080808080808080ULL);
        r12 = ((r12 - r16) + (r16 >> 7));
        r16 = (((r13 << 8) & 0x7f007f007f007f00ULL)
            | ((r13 & 0x7f007f007f007f00ULL) >> 8));
        r13 = ((r13 ^ 0x7f007f007f007f00ULL) + r16);
        r16 = (r13 & 0x8080808080808080ULL);
        r13 = ((r13 - r16) + (r16 >> 7));
        r16 = (((r14 << 8) & 0x7f007f007f007f00ULL)
            | ((r14 & 0x7f007f007f007f00ULL) >> 8));
        r14 = ((r14 ^ 0x7f007f007f007f00ULL) + r16);
        r16 = (r14 & 0x8080808080808080ULL);
        r14 = ((r14 - r16) + (r16 >> 7));
        r16 = (((r15 << 8) & 0x7f007f007f007f00ULL)
            | ((r15 & 0x7f007f007f007f00ULL) >> 8));
        r15 = ((r15 ^ 0x7f007f007f007f00ULL) + r16);
        r16 = (r15 & 0x8080808080808080ULL);
        r15 = ((r15 - r16) + (r16 >> 7));
        // Vector is now  r(i) for i = 
        //  0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
        // Butterfly: v[i], v[i+2] = v[i]+v[i+2], v[i]-v[i+2]
        r16 = (((r0 << 16) & 0x7f7f00007f7f0000ULL)
            | ((r0 & 0x7f7f00007f7f0000ULL) >> 16));
        r0 = ((r0 ^ 0x7f7f00007f7f0000ULL) + r16);
        r16 = (r0 & 0x8080808080808080ULL);
        r0 = ((r0 - r16) + (r16 >> 7));
        r16 = (((r1 << 16) & 0x7f7f00007f7f0000ULL)
            | ((r1 & 0x7f7f00007f7f0000ULL) >> 16));
        r1 = ((r1 ^ 0x7f7f00007f7f0000ULL) + r16);
        r16 = (r1 & 0x8080808080808080ULL);
        r1 = ((r1 - r16) + (r16 >> 7));
        r16 = (((r2 << 16) & 0x7f7f00007f7f0000ULL)
            | ((r2 & 0x7f7f00007f7f0000ULL) >> 16));
        r2 = ((r2 ^ 0x7f7f00007f7f0000ULL) + r16);
        r16 = (r2 & 0x8080808080808080ULL);
        r2 = ((r2 - r16) + (r16 >> 7));
        r16 = (((r3 << 16) & 0x7f7f00007f7f0000ULL)
            | ((r3 & 0x7f7f00007f7f0000ULL) >> 16));
        r3 = ((r3 ^ 0x7f7f00007f7f0000ULL) + r16);
        r16 = (r3 & 0x8080808080808080ULL);
        r3 = ((r3 - r16) + (r16 >> 7));
        r16 = (((r4 << 16) & 0x7f7f00007f7f0000ULL)
            | ((r4 & 0x7f7f00007f7f0000ULL) >> 16));
        r4 = ((r4 ^ 0x7f7f00007f7f0000ULL) + r16);
        r16 = (r4 & 0x8080808080808080ULL);
        r4 = ((r4 - r16) + (r16 >> 7));
        r16 = (((r5 << 16) & 0x7f7f00007f7f0000ULL)
            | ((r5 & 0x7f7f00007f7f0000ULL) >> 16));
        r5 = ((r5 ^ 0x7f7f00007f7f0000ULL) + r16);
        r16 = (r5 & 0x8080808080808080ULL);
        r5 = ((r5 - r16) + (r16 >> 7));
        r16 = (((r6 << 16) & 0x7f7f00007f7f0000ULL)
            | ((r6 & 0x7f7f00007f7f0000ULL) >> 16));
        r6 = ((r6 ^ 0x7f7f00007f7f0000ULL) + r16);
        r16 = (r6 & 0x8080808080808080ULL);
        r6 = ((r6 - r16) + (r16 >> 7));
        r16 = (((r7 << 16) & 0x7f7f00007f7f0000ULL)
            | ((r7 & 0x7f7f00007f7f0000ULL) >> 16));
        r7 = ((r7 ^ 0x7f7f00007f7f0000ULL) + r16);
        r16 = (r7 & 0x8080808080808080ULL);
        r7 = ((r7 - r16) + (r16 >> 7));
        r16 = (((r8 << 16) & 0x7f7f00007f7f0000ULL)
            | ((r8 & 0x7f7f00007f7f0000ULL) >> 16));
        r8 = ((r8 ^ 0x7f7f00007f7f0000ULL) + r16);
        r16 = (r8 & 0x8080808080808080ULL);
        r8 = ((r8 - r16) + (r16 >> 7));
        r16 = (((r9 << 16) & 0x7f7f00007f7f0000ULL)
            | ((r9 & 0x7f7f00007f7f0000ULL) >> 16));
        r9 = ((r9 ^ 0x7f7f00007f7f0000ULL) + r16);
        r16 = (r9 & 0x8080808080808080ULL);
        r9 = ((r9 - r16) + (r16 >> 7));
        r16 = (((r10 << 16) & 0x7f7f00007f7f0000ULL)
            | ((r10 & 0x7f7f00007f7f0000ULL) >> 16));
        r10 = ((r10 ^ 0x7f7f00007f7f0000ULL) + r16);
        r16 = (r10 & 0x8080808080808080ULL);
        r10 = ((r10 - r16) + (r16 >> 7));
        r16 = (((r11 << 16) & 0x7f7f00007f7f0000ULL)
            | ((r11 & 0x7f7f00007f7f0000ULL) >> 16));
        r11 = ((r11 ^ 0x7f7f00007f7f0000ULL) + r16);
        r16 = (r11 & 0x8080808080808080ULL);
        r11 = ((r11 - r16) + (r16 >> 7));
        r16 = (((r12 << 16) & 0x7f7f00007f7f0000ULL)
            | ((r12 & 0x7f7f00007f7f0000ULL) >> 16));
        r12 = ((r12 ^ 0x7f7f00007f7f0000ULL) + r16);
        r16 = (r12 & 0x8080808080808080ULL);
        r12 = ((r12 - r16) + (r16 >> 7));
        r16 = (((r13 << 16) & 0x7f7f00007f7f0000ULL)
            | ((r13 & 0x7f7f00007f7f0000ULL) >> 16));
        r13 = ((r13 ^ 0x7f7f00007f7f0000ULL) + r16);
        r16 = (r13 & 0x8080808080808080ULL);
        r13 = ((r13 - r16) + (r16 >> 7));
        r16 = (((r14 << 16) & 0x7f7f00007f7f0000ULL)
            | ((r14 & 0x7f7f00007f7f0000ULL) >> 16));
        r14 = ((r14 ^ 0x7f7f00007f7f0000ULL) + r16);
        r16 = (r14 & 0x8080808080808080ULL);
        r14 = ((r14 - r16) + (r16 >> 7));
        r16 = (((r15 << 16) & 0x7f7f00007f7f0000ULL)
            | ((r15 & 0x7f7f00007f7f0000ULL) >> 16));
        r15 = ((r15 ^ 0x7f7f00007f7f0000ULL) + r16);
        r16 = (r15 & 0x8080808080808080ULL);
        r15 = ((r15 - r16) + (r16 >> 7));
        // Vector is now  r(i) for i = 
        //  0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
        // Butterfly: v[i], v[i+8] = v[i]+v[i+8], v[i]-v[i+8]
        r16 = (r0 + (r1 ^ 0x7f7f7f7f7f7f7f7fULL));
        r0 = (r0 + r1);
        r1 = (r16 & 0x8080808080808080ULL);
        r1 = ((r16 - r1) + (r1 >> 7));
        r16 = (r0 & 0x8080808080808080ULL);
        r0 = ((r0 - r16) + (r16 >> 7));
        r16 = (r2 + (r3 ^ 0x7f7f7f7f7f7f7f7fULL));
        r2 = (r2 + r3);
        r3 = (r16 & 0x8080808080808080ULL);
        r3 = ((r16 - r3) + (r3 >> 7));
        r16 = (r2 & 0x8080808080808080ULL);
        r2 = ((r2 - r16) + (r16 >> 7));
        r16 = (r4 + (r5 ^ 0x7f7f7f7f7f7f7f7fULL));
        r4 = (r4 + r5);
        r5 = (r16 & 0x8080808080808080ULL);
        r5 = ((r16 - r5) + (r5 >> 7));
        r16 = (r4 & 0x8080808080808080ULL);
        r4 = ((r4 - r16) + (r16 >> 7));
        r16 = (r6 + (r7 ^ 0x7f7f7f7f7f7f7f7fULL));
        r6 = (r6 + r7);
        r7 = (r16 & 0x8080808080808080ULL);
        r7 = ((r16 - r7) + (r7 >> 7));
        r16 = (r6 & 0x8080808080808080ULL);
        r6 = ((r6 - r16) + (r16 >> 7));
        r16 = (r8 + (r9 ^ 0x7f7f7f7f7f7f7f7fULL));
        r8 = (r8 + r9);
        r9 = (r16 & 0x8080808080808080ULL);
        r9 = ((r16 - r9) + (r9 >> 7));
        r16 = (r8 & 0x8080808080808080ULL);
        r8 = ((r8 - r16) + (r16 >> 7));
        r16 = (r10 + (r11 ^ 0x7f7f7f7f7f7f7f7fULL));
        r10 = (r10 + r11);
        r11 = (r16 & 0x8080808080808080ULL);
        r11 = ((r16 - r11) + (r11 >> 7));
        r16 = (r10 & 0x8080808080808080ULL);
        r10 = ((r10 - r16) + (r16 >> 7));
        r16 = (r12 + (r13 ^ 0x7f7f7f7f7f7f7f7fULL));
        r12 = (r12 + r13);
        r13 = (r16 & 0x8080808080808080ULL);
        r13 = ((r16 - r13) + (r13 >> 7));
        r16 = (r12 & 0x8080808080808080ULL);
        r12 = ((r12 - r16) + (r16 >> 7));
        r16 = (r14 + (r15 ^ 0x7f7f7f7f7f7f7f7fULL));
        r14 = (r14 + r15);
        r15 = (r16 & 0x8080808080808080ULL);
        r15 = ((r16 - r15) + (r15 >> 7));
        r16 = (r14 & 0x8080808080808080ULL);
        r14 = ((r14 - r16) + (r16 >> 7));
        // Vector is now  r(i) for i = 
        //  0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
        // Butterfly: v[i], v[i+16] = v[i]+v[i+16], v[i]-v[i+16]
        r16 = (r0 + (r2 ^ 0x7f7f7f7f7f7f7f7fULL));
        r0 = (r0 + r2);
        r2 = (r16 & 0x8080808080808080ULL);
        r2 = ((r16 - r2) + (r2 >> 7));
        r16 = (r0 & 0x8080808080808080ULL);
        r0 = ((r0 - r16) + (r16 >> 7));
        r16 = (r1 + (r3 ^ 0x7f7f7f7f7f7f7f7fULL));
        r1 = (r1 + r3);
        r3 = (r16 & 0x8080808080808080ULL);
        r3 = ((r16 - r3) + (r3 >> 7));
        r16 = (r1 & 0x8080808080808080ULL);
        r1 = ((r1 - r16) + (r16 >> 7));
        r16 = (r4 + (r6 ^ 0x7f7f7f7f7f7f7f7fULL));
        r4 = (r4 + r6);
        r6 = (r16 & 0x8080808080808080ULL);
        r6 = ((r16 - r6) + (r6 >> 7));
        r16 = (r4 & 0x8080808080808080ULL);
        r4 = ((r4 - r16) + (r16 >> 7));
        r16 = (r5 + (r7 ^ 0x7f7f7f7f7f7f7f7fULL));
        r5 = (r5 + r7);
        r7 = (r16 & 0x8080808080808080ULL);
        r7 = ((r16 - r7) + (r7 >> 7));
        r16 = (r5 & 0x8080808080808080ULL);
        r5 = ((r5 - r16) + (r16 >> 7));
        r16 = (r8 + (r10 ^ 0x7f7f7f7f7f7f7f7fULL));
        r8 = (r8 + r10);
        r10 = (r16 & 0x8080808080808080ULL);
        r10 = ((r16 - r10) + (r10 >> 7));
        r16 = (r8 & 0x8080808080808080ULL);
        r8 = ((r8 - r16) + (r16 >> 7));
        r16 = (r9 + (r11 ^ 0x7f7f7f7f7f7f7f7fULL));
        r9 = (r9 + r11);
        r11 = (r16 & 0x8080808080808080ULL);
        r11 = ((r16 - r11) + (r11 >> 7));
        r16 = (r9 & 0x8080808080808080ULL);
        r9 = ((r9 - r16) + (r16 >> 7));
        r16 = (r12 + (r14 ^ 0x7f7f7f7f7f7f7f7fULL));
        r12 = (r12 + r14);
        r14 = (r16 & 0x8080808080808080ULL);
        r14 = ((r16 - r14) + (r14 >> 7));
        r16 = (r12 & 0x8080808080808080ULL);
        r12 = ((r12 - r16) + (r16 >> 7));
        r16 = (r13 + (r15 ^ 0x7f7f7f7f7f7f7f7fULL));
        r13 = (r13 + r15);
        r15 = (r16 & 0x8080808080808080ULL);
        r15 = ((r16 - r15) + (r15 >> 7));
        r16 = (r13 & 0x8080808080808080ULL);
        r13 = ((r13 - r16) + (r16 >> 7));
        // Vector is now  r(i) for i = 
        //  0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
        // Butterfly: v[i], v[i+32] = v[i]+v[i+32], v[i]-v[i+32]
        r16 = (r0 + (r4 ^ 0x7f7f7f7f7f7f7f7fULL));
        r0 = (r0 + r4);
        r4 = (r16 & 0x8080808080808080ULL);
        r4 = ((r16 - r4) + (r4 >> 7));
        r16 = (r0 & 0x8080808080808080ULL);
        r0 = ((r0 - r16) + (r16 >> 7));
        r16 = (r1 + (r5 ^ 0x7f7f7f7f7f7f7f7fULL));
        r1 = (r1 + r5);
        r5 = (r16 & 0x8080808080808080ULL);
        r5 = ((r16 - r5) + (r5 >> 7));
        r16 = (r1 & 0x8080808080808080ULL);
        r1 = ((r1 - r16) + (r16 >> 7));
        r16 = (r2 + (r6 ^ 0x7f7f7f7f7f7f7f7fULL));
        r2 = (r2 + r6);
        r6 = (r16 & 0x8080808080808080ULL);
        r6 = ((r16 - r6) + (r6 >> 7));
        r16 = (r2 & 0x8080808080808080ULL);
        r2 = ((r2 - r16) + (r16 >> 7));
        r16 = (r3 + (r7 ^ 0x7f7f7f7f7f7f7f7fULL));
        r3 = (r3 + r7);
        r7 = (r16 & 0x8080808080808080ULL);
        r7 = ((r16 - r7) + (r7 >> 7));
        r16 = (r3 & 0x8080808080808080ULL);
        r3 = ((r3 - r16) + (r16 >> 7));
        r16 = (r8 + (r12 ^ 0x7f7f7f7f7f7f7f7fULL));
        r8 = (r8 + r12);
        r12 = (r16 & 0x8080808080808080ULL);
        r12 = ((r16 - r12) + (r12 >> 7));
        r16 = (r8 & 0x8080808080808080ULL);
        r8 = ((r8 - r16) + (r16 >> 7));
        r16 = (r9 + (r13 ^ 0x7f7f7f7f7f7f7f7fULL));
        r9 = (r9 + r13);
        r13 = (r16 & 0x8080808080808080ULL);
        r13 = ((r16 - r13) + (r13 >> 7));
        r16 = (r9 & 0x8080808080808080ULL);
        r9 = ((r9 - r16) + (r16 >> 7));
        r16 = (r10 + (r14 ^ 0x7f7f7f7f7f7f7f7fULL));
        r10 = (r10 + r14);
        r14 = (r16 & 0x8080808080808080ULL);
        r14 = ((r16 - r14) + (r14 >> 7));
        r16 = (r10 & 0x8080808080808080ULL);
        r10 = ((r10 - r16) + (r16 >> 7));
        r16 = (r11 + (r15 ^ 0x7f7f7f7f7f7f7f7fULL));
        r11 = (r11 + r15);
        r15 = (r16 & 0x8080808080808080ULL);
        r15 = ((r16 - r15) + (r15 >> 7));
        r16 = (r11 & 0x8080808080808080ULL);
        r11 = ((r11 - r16) + (r16 >> 7));
        // Vector is now  r(i) for i = 
        //  0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
        // Butterfly: v[i], v[i+64] = v[i]+v[i+64], v[i]-v[i+64]
        r16 = (r0 + (r8 ^ 0x7f7f7f7f7f7f7f7fULL));
        r0 = (r0 + r8);
        r8 = (r16 & 0x8080808080808080ULL);
        r8 = ((r16 - r8) + (r8 >> 7));
        r16 = (r0 & 0x8080808080808080ULL);
        r0 = ((r0 - r16) + (r16 >> 7));
        r16 = (r1 + (r9 ^ 0x7f7f7f7f7f7f7f7fULL));
        r1 = (r1 + r9);
        r9 = (r16 & 0x8080808080808080ULL);
        r9 = ((r16 - r9) + (r9 >> 7));
        r16 = (r1 & 0x8080808080808080ULL);
        r1 = ((r1 - r16) + (r16 >> 7));
        r16 = (r2 + (r10 ^ 0x7f7f7f7f7f7f7f7fULL));
        r2 = (r2 + r10);
        r10 = (r16 & 0x8080808080808080ULL);
        r10 = ((r16 - r10) + (r10 >> 7));
        r16 = (r2 & 0x8080808080808080ULL);
        r2 = ((r2 - r16) + (r16 >> 7));
        r16 = (r3 + (r11 ^ 0x7f7f7f7f7f7f7f7fULL));
        r3 = (r3 + r11);
        r11 = (r16 & 0x8080808080808080ULL);
        r11 = ((r16 - r11) + (r11 >> 7));
        r16 = (r3 & 0x8080808080808080ULL);
        r3 = ((r3 - r16) + (r16 >> 7));
        r16 = (r4 + (r12 ^ 0x7f7f7f7f7f7f7f7fULL));
        r4 = (r4 + r12);
        r12 = (r16 & 0x8080808080808080ULL);
        r12 = ((r16 - r12) + (r12 >> 7));
        r16 = (r4 & 0x8080808080808080ULL);
        r4 = ((r4 - r16) + (r16 >> 7));
        r16 = (r5 + (r13 ^ 0x7f7f7f7f7f7f7f7fULL));
        r5 = (r5 + r13);
        r13 = (r16 & 0x8080808080808080ULL);
        r13 = ((r16 - r13) + (r13 >> 7));
        r16 = (r5 & 0x8080808080808080ULL);
        r5 = ((r5 - r16) + (r16 >> 7));
        r16 = (r6 + (r14 ^ 0x7f7f7f7f7f7f7f7fULL));
        r6 = (r6 + r14);
        r14 = (r16 & 0x8080808080808080ULL);
        r14 = ((r16 - r14) + (r14 >> 7));
        r16 = (r6 & 0x8080808080808080ULL);
        r6 = ((r6 - r16) + (r16 >> 7));
        r16 = (r7 + (r15 ^ 0x7f7f7f7f7f7f7f7fULL));
        r7 = (r7 + r15);
        r15 = (r16 & 0x8080808080808080ULL);
        r15 = ((r16 - r15) + (r15 >> 7));
        r16 = (r7 & 0x8080808080808080ULL);
        r7 = ((r7 - r16) + (r16 >> 7));
        // Vector is now  r(i) for i = 
        //  0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
        r0 = (((r0 & 0x707070707070707ULL) << 4)
            | ((r0 & 0x7878787878787878ULL) >> 3));
        v_out[0] = r0 ^  p_mask[6];
        r1 = (((r1 & 0x707070707070707ULL) << 4)
            | ((r1 & 0x7878787878787878ULL) >> 3));
        v_out[4] = r1 ^  p_mask[6];
        r2 = (((r2 & 0x707070707070707ULL) << 4)
            | ((r2 & 0x7878787878787878ULL) >> 3));
        v_out[8] = r2 ^  p_mask[6];
        r3 = (((r3 & 0x707070707070707ULL) << 4)
            | ((r3 & 0x7878787878787878ULL) >> 3));
        v_out[12] = r3 ^  p_mask[4];
        r4 = (((r4 & 0x707070707070707ULL) << 4)
            | ((r4 & 0x7878787878787878ULL) >> 3));
        v_out[16] = r4 ^  p_mask[6];
        r5 = (((r5 & 0x707070707070707ULL) << 4)
            | ((r5 & 0x7878787878787878ULL) >> 3));
        v_out[20] = r5 ^  p_mask[4];
        r6 = (((r6 & 0x707070707070707ULL) << 4)
            | ((r6 & 0x7878787878787878ULL) >> 3));
        v_out[24] = r6 ^  p_mask[4];
        r7 = (((r7 & 0x707070707070707ULL) << 4)
            | ((r7 & 0x7878787878787878ULL) >> 3));
        v_out[28] = r7 ^  p_mask[4];
        r8 = (((r8 & 0x707070707070707ULL) << 4)
            | ((r8 & 0x7878787878787878ULL) >> 3));
        v_out[32] = r8 ^  p_mask[6];
        r9 = (((r9 & 0x707070707070707ULL) << 4)
            | ((r9 & 0x7878787878787878ULL) >> 3));
        v_out[36] = r9 ^  p_mask[4];
        r10 = (((r10 & 0x707070707070707ULL) << 4)
            | ((r10 & 0x7878787878787878ULL) >> 3));
        v_out[40] = r10 ^  p_mask[4];
        r11 = (((r11 & 0x707070707070707ULL) << 4)
            | ((r11 & 0x7878787878787878ULL) >> 3));
        v_out[44] = r11 ^  p_mask[4];
        r12 = (((r12 & 0x707070707070707ULL) << 4)
            | ((r12 & 0x7878787878787878ULL) >> 3));
        v_out[48] = r12 ^  p_mask[4];
        r13 = (((r13 & 0x707070707070707ULL) << 4)
            | ((r13 & 0x7878787878787878ULL) >> 3));
        v_out[52] = r13 ^  p_mask[4];
        r14 = (((r14 & 0x707070707070707ULL) << 4)
            | ((r14 & 0x7878787878787878ULL) >> 3));
        v_out[56] = r14 ^  p_mask[4];
        r15 = (((r15 & 0x707070707070707ULL) << 4)
            | ((r15 & 0x7878787878787878ULL) >> 3));
        v_out[60] = r15 ^  p_mask[6];
        v_in++;
        v_out++;
        }
        v_out[0] = 0;
        v_out[4] = 0;
        v_out[8] = 0;
        v_out[12] = 0;
        v_out[16] = 0;
        v_out[20] = 0;
        v_out[24] = 0;
        v_out[28] = 0;
        v_out[32] = 0;
        v_out[36] = 0;
        v_out[40] = 0;
        v_out[44] = 0;
        v_out[48] = 0;
        v_out[52] = 0;
        v_out[56] = 0;
        v_out[60] = 0;
        v_in += 61;
        v_out += 61;
        }
        // End of automatically generated matrix operation.
 
    }
}  


static void mm_op127_xi_a(uint_mmv_t *v_in,  uint32_t exp1, uint_mmv_t *v_out)
{
    uint_fast32_t i;
    uint_mmv_t e_mask =  -((uint_mmv_t)exp1 & 0x1ULL);
    for (i = 0; i < 6; ++i) {
        // %%MUL_MATRIX_XI16 v_in, e_mask, v_out

        // This is an automatically generated matrix operation, do not change!
        {
        uint_mmv_t r0, r1, r2, r3, r4;

        uint_fast32_t i;
        // TODO: write comment!!!
        // 
        for (i = 0; i < 3; ++i) {
        e_mask = ~(e_mask);
        r0 = v_in[0] ^  (0x7f7f7f007f7f7f00ULL & e_mask);
        r4 = ((r0 ^ (r0 >> 8)) & 0x7f0000007f00ULL);
        r0 ^= (r4 | (r4 << 8));
        r1 = v_in[8] ^  (0x7f0000007fULL & e_mask);
        r4 = ((r1 ^ (r1 >> 8)) & 0x7f0000007f00ULL);
        r1 ^= (r4 | (r4 << 8));
        r2 = v_in[4] ^  (0x7f0000007fULL & e_mask);
        r4 = ((r2 ^ (r2 >> 8)) & 0x7f0000007f00ULL);
        r2 ^= (r4 | (r4 << 8));
        r3 = v_in[12] ^  (0x7f0000007fULL & e_mask);
        r4 = ((r3 ^ (r3 >> 8)) & 0x7f0000007f00ULL);
        r3 ^= (r4 | (r4 << 8));
        // Butterfly: v[i], v[i+1] = v[i]+v[i+1], v[i]-v[i+1]
        r4 = (((r0 << 8) & 0x7f007f007f007f00ULL)
            | ((r0 & 0x7f007f007f007f00ULL) >> 8));
        r0 = ((r0 ^ 0x7f007f007f007f00ULL) + r4);
        r4 = (r0 & 0x8080808080808080ULL);
        r0 = ((r0 - r4) + (r4 >> 7));
        r4 = (((r1 << 8) & 0x7f007f007f007f00ULL)
            | ((r1 & 0x7f007f007f007f00ULL) >> 8));
        r1 = ((r1 ^ 0x7f007f007f007f00ULL) + r4);
        r4 = (r1 & 0x8080808080808080ULL);
        r1 = ((r1 - r4) + (r4 >> 7));
        r4 = (((r2 << 8) & 0x7f007f007f007f00ULL)
            | ((r2 & 0x7f007f007f007f00ULL) >> 8));
        r2 = ((r2 ^ 0x7f007f007f007f00ULL) + r4);
        r4 = (r2 & 0x8080808080808080ULL);
        r2 = ((r2 - r4) + (r4 >> 7));
        r4 = (((r3 << 8) & 0x7f007f007f007f00ULL)
            | ((r3 & 0x7f007f007f007f00ULL) >> 8));
        r3 = ((r3 ^ 0x7f007f007f007f00ULL) + r4);
        r4 = (r3 & 0x8080808080808080ULL);
        r3 = ((r3 - r4) + (r4 >> 7));
        // Vector is now  r(i) for i = 0,1,2,3
        // Butterfly: v[i], v[i+2] = v[i]+v[i+2], v[i]-v[i+2]
        r4 = (((r0 << 16) & 0x7f7f00007f7f0000ULL)
            | ((r0 & 0x7f7f00007f7f0000ULL) >> 16));
        r0 = ((r0 ^ 0x7f7f00007f7f0000ULL) + r4);
        r4 = (r0 & 0x8080808080808080ULL);
        r0 = ((r0 - r4) + (r4 >> 7));
        r4 = (((r1 << 16) & 0x7f7f00007f7f0000ULL)
            | ((r1 & 0x7f7f00007f7f0000ULL) >> 16));
        r1 = ((r1 ^ 0x7f7f00007f7f0000ULL) + r4);
        r4 = (r1 & 0x8080808080808080ULL);
        r1 = ((r1 - r4) + (r4 >> 7));
        r4 = (((r2 << 16) & 0x7f7f00007f7f0000ULL)
            | ((r2 & 0x7f7f00007f7f0000ULL) >> 16));
        r2 = ((r2 ^ 0x7f7f00007f7f0000ULL) + r4);
        r4 = (r2 & 0x8080808080808080ULL);
        r2 = ((r2 - r4) + (r4 >> 7));
        r4 = (((r3 << 16) & 0x7f7f00007f7f0000ULL)
            | ((r3 & 0x7f7f00007f7f0000ULL) >> 16));
        r3 = ((r3 ^ 0x7f7f00007f7f0000ULL) + r4);
        r4 = (r3 & 0x8080808080808080ULL);
        r3 = ((r3 - r4) + (r4 >> 7));
        // Vector is now  r(i) for i = 0,1,2,3
        // Butterfly: v[i], v[i+8] = v[i]+v[i+8], v[i]-v[i+8]
        r4 = (r0 + (r1 ^ 0x7f7f7f7f7f7f7f7fULL));
        r0 = (r0 + r1);
        r1 = (r4 & 0x8080808080808080ULL);
        r1 = ((r4 - r1) + (r1 >> 7));
        r4 = (r0 & 0x8080808080808080ULL);
        r0 = ((r0 - r4) + (r4 >> 7));
        r4 = (r2 + (r3 ^ 0x7f7f7f7f7f7f7f7fULL));
        r2 = (r2 + r3);
        r3 = (r4 & 0x8080808080808080ULL);
        r3 = ((r4 - r3) + (r3 >> 7));
        r4 = (r2 & 0x8080808080808080ULL);
        r2 = ((r2 - r4) + (r4 >> 7));
        // Vector is now  r(i) for i = 0,1,2,3
        // Butterfly: v[i], v[i+16] = v[i]+v[i+16], v[i]-v[i+16]
        r4 = (r0 + (r2 ^ 0x7f7f7f7f7f7f7f7fULL));
        r0 = (r0 + r2);
        r2 = (r4 & 0x8080808080808080ULL);
        r2 = ((r4 - r2) + (r2 >> 7));
        r4 = (r0 & 0x8080808080808080ULL);
        r0 = ((r0 - r4) + (r4 >> 7));
        r4 = (r1 + (r3 ^ 0x7f7f7f7f7f7f7f7fULL));
        r1 = (r1 + r3);
        r3 = (r4 & 0x8080808080808080ULL);
        r3 = ((r4 - r3) + (r3 >> 7));
        r4 = (r1 & 0x8080808080808080ULL);
        r1 = ((r1 - r4) + (r4 >> 7));
        // Vector is now  r(i) for i = 0,1,2,3
        e_mask = ~(e_mask);
        r0 = (((r0 & 0x303030303030303ULL) << 5)
            | ((r0 & 0x7c7c7c7c7c7c7c7cULL) >> 2));
        v_out[0] = r0 ^ (e_mask & 0x7f7f7f007f7f7f00ULL);
        r1 = (((r1 & 0x303030303030303ULL) << 5)
            | ((r1 & 0x7c7c7c7c7c7c7c7cULL) >> 2));
        v_out[4] = r1 ^ (e_mask & 0x7f0000007fULL);
        r2 = (((r2 & 0x303030303030303ULL) << 5)
            | ((r2 & 0x7c7c7c7c7c7c7c7cULL) >> 2));
        v_out[8] = r2 ^ (e_mask & 0x7f0000007fULL);
        r3 = (((r3 & 0x303030303030303ULL) << 5)
            | ((r3 & 0x7c7c7c7c7c7c7c7cULL) >> 2));
        v_out[12] = r3 ^ (e_mask & 0x7f0000007fULL);
        // 78 lines of code, 194 operations
        v_in++;
        v_out++;
        }
        v_out[0] = 0;
        v_out[4] = 0;
        v_out[8] = 0;
        v_out[12] = 0;
        v_in += 13;
        v_out += 13;
        }
        // End of automatically generated matrix operation.
 
    }
}  


// %%EXPORT p
void mm_op127_xi(uint_mmv_t *v_in,  uint32_t exp, uint_mmv_t *v_out)
{
    uint_mmv_t i, exp1;
 
    exp %= 3;
    if (exp == 0) {
        for (i = 0; i < 30936; ++i) v_out[i] = v_in[i];
        return;
    }
    exp1 =  (uint_mmv_t)exp - 0x1ULL;

    // Do monomial part, i.e. tags B, C, T, X
    mm_op127_xi_mon(v_in, exp1, v_out);

    // Do tag A
    mm_op127_xi_a(v_in, exp1, v_out); 

    // Do tags X, Y
    for (i = 0; i < 4; ++i) {
        uint_mmv_t *p_src = v_in + MM_OP127_OFS_Z + (i << HALF_YZ_SHIFT);
        mm_op127_xi_yz(p_src, exp1, v_out + TAB127_XI64_OFFSET[exp1][i]);
    }
    
}


