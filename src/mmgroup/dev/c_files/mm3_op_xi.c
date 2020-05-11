/////////////////////////////////////////////////////////////////////////////
// This C file has been created automatically. Do not edit!!!
/////////////////////////////////////////////////////////////////////////////

// %%COMMENT
// TODO: Yet to be documented!!!


#include <string.h>
#include "mat24_functions.h"
#include "mm_op3.h"   


static void mm_op3_xi_mon(uint_mmv_t *v_in,  uint32_t exp1, uint_mmv_t *v_out)
{
    uint_mmv_t *p_src, *p_dest;
    uint_fast32_t i, j;
    uint_fast32_t diff = exp1 ? 1024 : 0;
    uint8_t b[2496], *p_b;
    mm_sub_table_xi_type *p_tables = mm_sub_table_xi[exp1];
    uint16_t *p_perm;
    uint32_t *p_sign;



    ///////////////////////////////////////////////////////////////
    // Map tag BC to tag BC.
    ///////////////////////////////////////////////////////////////
    p_src = v_in + 24;
    p_dest = v_out + 24;
    p_sign = p_tables[0].p_sign;
    p_perm = p_tables[0].p_perm;

    for (i = 0; i < 1; ++i) {
        p_b = b;
        for (j = 0; j < 78; ++j) {
           // %%OP_XI_LOAD p_src, p_b, {i_src.SHAPE[2]}
           uint_mmv_t r0, r1, r2, r3;
           r0 =  (p_src)[0];
           r1 = (r0 >> 2) & 0x303030303030303ULL;
           r2 = (r0 >> 4) & 0x303030303030303ULL;
           r3 = (r0 >> 6) & 0x303030303030303ULL;
           r0 &= 0x303030303030303ULL;
           (p_b)[0] = r0 >> 0;
           (p_b)[1] = r1 >> 0;
           (p_b)[2] = r2 >> 0;
           (p_b)[3] = r3 >> 0;
           (p_b)[4] = r0 >> 8;
           (p_b)[5] = r1 >> 8;
           (p_b)[6] = r2 >> 8;
           (p_b)[7] = r3 >> 8;
           (p_b)[8] = r0 >> 16;
           (p_b)[9] = r1 >> 16;
           (p_b)[10] = r2 >> 16;
           (p_b)[11] = r3 >> 16;
           (p_b)[12] = r0 >> 24;
           (p_b)[13] = r1 >> 24;
           (p_b)[14] = r2 >> 24;
           (p_b)[15] = r3 >> 24;
           (p_b)[16] = r0 >> 32;
           (p_b)[17] = r1 >> 32;
           (p_b)[18] = r2 >> 32;
           (p_b)[19] = r3 >> 32;
           (p_b)[20] = r0 >> 40;
           (p_b)[21] = r1 >> 40;
           (p_b)[22] = r2 >> 40;
           (p_b)[23] = r3 >> 40;
           (p_b)[24] = r0 >> 48;
           (p_b)[25] = r1 >> 48;
           (p_b)[26] = r2 >> 48;
           (p_b)[27] = r3 >> 48;
           (p_b)[28] = r0 >> 56;
           (p_b)[29] = r1 >> 56;
           (p_b)[30] = r2 >> 56;
           (p_b)[31] = r3 >> 56;
           p_src += 1;
           p_b += 32;
        }

        for (j = 0; j < 78; ++j) {
           // %%OP_XI_STORE b, p_perm, p_sign, p_dest, {i_dest.SHAPE[2]}
           uint_mmv_t r0, r1;
           r0 = ((uint_mmv_t)(b[p_perm[0]]) << 0)
             + ((uint_mmv_t)(b[p_perm[1]]) << 2)
             + ((uint_mmv_t)(b[p_perm[2]]) << 4)
             + ((uint_mmv_t)(b[p_perm[3]]) << 6)
             + ((uint_mmv_t)(b[p_perm[4]]) << 8)
             + ((uint_mmv_t)(b[p_perm[5]]) << 10)
             + ((uint_mmv_t)(b[p_perm[6]]) << 12)
             + ((uint_mmv_t)(b[p_perm[7]]) << 14)
             + ((uint_mmv_t)(b[p_perm[8]]) << 16)
             + ((uint_mmv_t)(b[p_perm[9]]) << 18)
             + ((uint_mmv_t)(b[p_perm[10]]) << 20)
             + ((uint_mmv_t)(b[p_perm[11]]) << 22)
             + ((uint_mmv_t)(b[p_perm[12]]) << 24)
             + ((uint_mmv_t)(b[p_perm[13]]) << 26)
             + ((uint_mmv_t)(b[p_perm[14]]) << 28)
             + ((uint_mmv_t)(b[p_perm[15]]) << 30)
             + ((uint_mmv_t)(b[p_perm[16]]) << 32)
             + ((uint_mmv_t)(b[p_perm[17]]) << 34)
             + ((uint_mmv_t)(b[p_perm[18]]) << 36)
             + ((uint_mmv_t)(b[p_perm[19]]) << 38)
             + ((uint_mmv_t)(b[p_perm[20]]) << 40)
             + ((uint_mmv_t)(b[p_perm[21]]) << 42)
             + ((uint_mmv_t)(b[p_perm[22]]) << 44)
             + ((uint_mmv_t)(b[p_perm[23]]) << 46)
             + ((uint_mmv_t)(b[p_perm[24]]) << 48)
             + ((uint_mmv_t)(b[p_perm[25]]) << 50)
             + ((uint_mmv_t)(b[p_perm[26]]) << 52)
             + ((uint_mmv_t)(b[p_perm[27]]) << 54)
             + ((uint_mmv_t)(b[p_perm[28]]) << 56)
             + ((uint_mmv_t)(b[p_perm[29]]) << 58)
             + ((uint_mmv_t)(b[p_perm[30]]) << 60)
             + ((uint_mmv_t)(b[p_perm[31]]) << 62);
           r1 = (p_sign)[0] >> 0;
           // Spread bits 0,...,31 of r1 to the (2-bit long) fields
           // of r1. A field of r1 is set to 0x3 if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xffffULL) 
               +  ((r1 & 0xffff0000ULL) << 16);
           r1 = (r1 & 0xff000000ffULL) 
               +  ((r1 & 0xff000000ff00ULL) << 8);
           r1 = (r1 & 0xf000f000f000fULL) 
               +  ((r1 & 0xf000f000f000f0ULL) << 4);
           r1 = (r1 & 0x303030303030303ULL) 
               +  ((r1 & 0xc0c0c0c0c0c0c0cULL) << 2);
           r1 = (r1 & 0x1111111111111111ULL) 
               +  ((r1 & 0x2222222222222222ULL) << 1);
           r1 = (((r1) << 2) - (r1));
           // Bit spreading done.
           (p_dest)[0] = r0 ^ r1;
           p_dest += 1;
           p_perm += 32;
           p_sign += 1;
        }
    }

    ///////////////////////////////////////////////////////////////
    // Map tag T0 to tag T0.
    ///////////////////////////////////////////////////////////////
    p_src = v_in + 102;
    p_dest = v_out + 102;
    p_sign = p_tables[1].p_sign;
    p_perm = p_tables[1].p_perm;

    for (i = 0; i < 45; ++i) {
        p_b = b;
        for (j = 0; j < 16; ++j) {
           // %%OP_XI_LOAD p_src, p_b, {i_src.SHAPE[2]}
           uint_mmv_t r0, r1, r2, r3;
           r0 =  (p_src)[0];
           r1 = (r0 >> 2) & 0x303030303030303ULL;
           r2 = (r0 >> 4) & 0x303030303030303ULL;
           r3 = (r0 >> 6) & 0x303030303030303ULL;
           r0 &= 0x303030303030303ULL;
           (p_b)[0] = r0 >> 0;
           (p_b)[1] = r1 >> 0;
           (p_b)[2] = r2 >> 0;
           (p_b)[3] = r3 >> 0;
           (p_b)[4] = r0 >> 8;
           (p_b)[5] = r1 >> 8;
           (p_b)[6] = r2 >> 8;
           (p_b)[7] = r3 >> 8;
           (p_b)[8] = r0 >> 16;
           (p_b)[9] = r1 >> 16;
           (p_b)[10] = r2 >> 16;
           (p_b)[11] = r3 >> 16;
           (p_b)[12] = r0 >> 24;
           (p_b)[13] = r1 >> 24;
           (p_b)[14] = r2 >> 24;
           (p_b)[15] = r3 >> 24;
           (p_b)[16] = r0 >> 32;
           (p_b)[17] = r1 >> 32;
           (p_b)[18] = r2 >> 32;
           (p_b)[19] = r3 >> 32;
           (p_b)[20] = r0 >> 40;
           (p_b)[21] = r1 >> 40;
           (p_b)[22] = r2 >> 40;
           (p_b)[23] = r3 >> 40;
           (p_b)[24] = r0 >> 48;
           (p_b)[25] = r1 >> 48;
           (p_b)[26] = r2 >> 48;
           (p_b)[27] = r3 >> 48;
           (p_b)[28] = r0 >> 56;
           (p_b)[29] = r1 >> 56;
           (p_b)[30] = r2 >> 56;
           (p_b)[31] = r3 >> 56;
           p_src += 1;
           p_b += 32;
        }

        for (j = 0; j < 16; ++j) {
           // %%OP_XI_STORE b, p_perm, p_sign, p_dest, {i_dest.SHAPE[2]}
           uint_mmv_t r0, r1;
           r0 = ((uint_mmv_t)(b[p_perm[0]]) << 0)
             + ((uint_mmv_t)(b[p_perm[1]]) << 2)
             + ((uint_mmv_t)(b[p_perm[2]]) << 4)
             + ((uint_mmv_t)(b[p_perm[3]]) << 6)
             + ((uint_mmv_t)(b[p_perm[4]]) << 8)
             + ((uint_mmv_t)(b[p_perm[5]]) << 10)
             + ((uint_mmv_t)(b[p_perm[6]]) << 12)
             + ((uint_mmv_t)(b[p_perm[7]]) << 14)
             + ((uint_mmv_t)(b[p_perm[8]]) << 16)
             + ((uint_mmv_t)(b[p_perm[9]]) << 18)
             + ((uint_mmv_t)(b[p_perm[10]]) << 20)
             + ((uint_mmv_t)(b[p_perm[11]]) << 22)
             + ((uint_mmv_t)(b[p_perm[12]]) << 24)
             + ((uint_mmv_t)(b[p_perm[13]]) << 26)
             + ((uint_mmv_t)(b[p_perm[14]]) << 28)
             + ((uint_mmv_t)(b[p_perm[15]]) << 30)
             + ((uint_mmv_t)(b[p_perm[16]]) << 32)
             + ((uint_mmv_t)(b[p_perm[17]]) << 34)
             + ((uint_mmv_t)(b[p_perm[18]]) << 36)
             + ((uint_mmv_t)(b[p_perm[19]]) << 38)
             + ((uint_mmv_t)(b[p_perm[20]]) << 40)
             + ((uint_mmv_t)(b[p_perm[21]]) << 42)
             + ((uint_mmv_t)(b[p_perm[22]]) << 44)
             + ((uint_mmv_t)(b[p_perm[23]]) << 46)
             + ((uint_mmv_t)(b[p_perm[24]]) << 48)
             + ((uint_mmv_t)(b[p_perm[25]]) << 50)
             + ((uint_mmv_t)(b[p_perm[26]]) << 52)
             + ((uint_mmv_t)(b[p_perm[27]]) << 54)
             + ((uint_mmv_t)(b[p_perm[28]]) << 56)
             + ((uint_mmv_t)(b[p_perm[29]]) << 58)
             + ((uint_mmv_t)(b[p_perm[30]]) << 60)
             + ((uint_mmv_t)(b[p_perm[31]]) << 62);
           r1 = (p_sign)[0] >> 0;
           // Spread bits 0,...,31 of r1 to the (2-bit long) fields
           // of r1. A field of r1 is set to 0x3 if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xffffULL) 
               +  ((r1 & 0xffff0000ULL) << 16);
           r1 = (r1 & 0xff000000ffULL) 
               +  ((r1 & 0xff000000ff00ULL) << 8);
           r1 = (r1 & 0xf000f000f000fULL) 
               +  ((r1 & 0xf000f000f000f0ULL) << 4);
           r1 = (r1 & 0x303030303030303ULL) 
               +  ((r1 & 0xc0c0c0c0c0c0c0cULL) << 2);
           r1 = (r1 & 0x1111111111111111ULL) 
               +  ((r1 & 0x2222222222222222ULL) << 1);
           r1 = (((r1) << 2) - (r1));
           // Bit spreading done.
           (p_dest)[0] = r0 ^ r1;
           p_dest += 1;
           p_perm += 32;
           p_sign += 1;
        }
    }

    ///////////////////////////////////////////////////////////////
    // Map tag T1 to tag X0 if e = 1
    // Map tag T1 to tag X1 if e = 2
    ///////////////////////////////////////////////////////////////
    p_src = v_in + 822;
    p_dest = v_out + 1590;
    p_dest += diff;
    p_sign = p_tables[2].p_sign;
    p_perm = p_tables[2].p_perm;

    for (i = 0; i < 64; ++i) {
        p_b = b;
        for (j = 0; j < 12; ++j) {
           // %%OP_XI_LOAD p_src, p_b, {i_src.SHAPE[2]}
           uint_mmv_t r0, r1, r2, r3;
           r0 =  (p_src)[0];
           r1 = (r0 >> 2) & 0x303030303030303ULL;
           r2 = (r0 >> 4) & 0x303030303030303ULL;
           r3 = (r0 >> 6) & 0x303030303030303ULL;
           r0 &= 0x303030303030303ULL;
           (p_b)[0] = r0 >> 0;
           (p_b)[1] = r1 >> 0;
           (p_b)[2] = r2 >> 0;
           (p_b)[3] = r3 >> 0;
           (p_b)[4] = r0 >> 8;
           (p_b)[5] = r1 >> 8;
           (p_b)[6] = r2 >> 8;
           (p_b)[7] = r3 >> 8;
           (p_b)[8] = r0 >> 16;
           (p_b)[9] = r1 >> 16;
           (p_b)[10] = r2 >> 16;
           (p_b)[11] = r3 >> 16;
           (p_b)[12] = r0 >> 24;
           (p_b)[13] = r1 >> 24;
           (p_b)[14] = r2 >> 24;
           (p_b)[15] = r3 >> 24;
           (p_b)[16] = r0 >> 32;
           (p_b)[17] = r1 >> 32;
           (p_b)[18] = r2 >> 32;
           (p_b)[19] = r3 >> 32;
           (p_b)[20] = r0 >> 40;
           (p_b)[21] = r1 >> 40;
           (p_b)[22] = r2 >> 40;
           (p_b)[23] = r3 >> 40;
           (p_b)[24] = r0 >> 48;
           (p_b)[25] = r1 >> 48;
           (p_b)[26] = r2 >> 48;
           (p_b)[27] = r3 >> 48;
           (p_b)[28] = r0 >> 56;
           (p_b)[29] = r1 >> 56;
           (p_b)[30] = r2 >> 56;
           (p_b)[31] = r3 >> 56;
           p_src += 1;
           p_b += 32;
        }

        for (j = 0; j < 16; ++j) {
           // %%OP_XI_STORE b, p_perm, p_sign, p_dest, {i_dest.SHAPE[2]}
           uint_mmv_t r0, r1;
           r0 = ((uint_mmv_t)(b[p_perm[0]]) << 0)
             + ((uint_mmv_t)(b[p_perm[1]]) << 2)
             + ((uint_mmv_t)(b[p_perm[2]]) << 4)
             + ((uint_mmv_t)(b[p_perm[3]]) << 6)
             + ((uint_mmv_t)(b[p_perm[4]]) << 8)
             + ((uint_mmv_t)(b[p_perm[5]]) << 10)
             + ((uint_mmv_t)(b[p_perm[6]]) << 12)
             + ((uint_mmv_t)(b[p_perm[7]]) << 14)
             + ((uint_mmv_t)(b[p_perm[8]]) << 16)
             + ((uint_mmv_t)(b[p_perm[9]]) << 18)
             + ((uint_mmv_t)(b[p_perm[10]]) << 20)
             + ((uint_mmv_t)(b[p_perm[11]]) << 22)
             + ((uint_mmv_t)(b[p_perm[12]]) << 24)
             + ((uint_mmv_t)(b[p_perm[13]]) << 26)
             + ((uint_mmv_t)(b[p_perm[14]]) << 28)
             + ((uint_mmv_t)(b[p_perm[15]]) << 30)
             + ((uint_mmv_t)(b[p_perm[16]]) << 32)
             + ((uint_mmv_t)(b[p_perm[17]]) << 34)
             + ((uint_mmv_t)(b[p_perm[18]]) << 36)
             + ((uint_mmv_t)(b[p_perm[19]]) << 38)
             + ((uint_mmv_t)(b[p_perm[20]]) << 40)
             + ((uint_mmv_t)(b[p_perm[21]]) << 42)
             + ((uint_mmv_t)(b[p_perm[22]]) << 44)
             + ((uint_mmv_t)(b[p_perm[23]]) << 46);
           r1 = (p_sign)[0] >> 0;
           // Spread bits 0,...,31 of r1 to the (2-bit long) fields
           // of r1. A field of r1 is set to 0x3 if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xffffULL) 
               +  ((r1 & 0xffff0000ULL) << 16);
           r1 = (r1 & 0xff000000ffULL) 
               +  ((r1 & 0xff000000ff00ULL) << 8);
           r1 = (r1 & 0xf000f000f000fULL) 
               +  ((r1 & 0xf000f000f000f0ULL) << 4);
           r1 = (r1 & 0x303030303030303ULL) 
               +  ((r1 & 0xc0c0c0c0c0c0c0cULL) << 2);
           r1 = (r1 & 0x1111111111111111ULL) 
               +  ((r1 & 0x2222222222222222ULL) << 1);
           r1 = (((r1) << 2) - (r1));
           // Bit spreading done.
           (p_dest)[0] = r0 ^ r1;
           p_dest += 1;
           p_perm += 24;
           p_sign += 1;
        }
    }

    ///////////////////////////////////////////////////////////////
    // Map tag X0 to tag X1 if e = 1
    // Map tag X1 to tag X0 if e = 2
    ///////////////////////////////////////////////////////////////
    p_src = v_in + 1590;
    p_src += diff;
    p_dest = v_out + 2614;
    p_dest -= diff;
    p_sign = p_tables[3].p_sign;
    p_perm = p_tables[3].p_perm;

    for (i = 0; i < 64; ++i) {
        p_b = b;
        for (j = 0; j < 16; ++j) {
           // %%OP_XI_LOAD p_src, p_b, {i_src.SHAPE[2]}
           uint_mmv_t r0, r1, r2, r3;
           r0 =  (p_src)[0];
           r1 = (r0 >> 2) & 0x303030303030303ULL;
           r2 = (r0 >> 4) & 0x303030303030303ULL;
           r3 = (r0 >> 6) & 0x303030303030303ULL;
           r0 &= 0x303030303030303ULL;
           (p_b)[0] = r0 >> 0;
           (p_b)[1] = r1 >> 0;
           (p_b)[2] = r2 >> 0;
           (p_b)[3] = r3 >> 0;
           (p_b)[4] = r0 >> 8;
           (p_b)[5] = r1 >> 8;
           (p_b)[6] = r2 >> 8;
           (p_b)[7] = r3 >> 8;
           (p_b)[8] = r0 >> 16;
           (p_b)[9] = r1 >> 16;
           (p_b)[10] = r2 >> 16;
           (p_b)[11] = r3 >> 16;
           (p_b)[12] = r0 >> 24;
           (p_b)[13] = r1 >> 24;
           (p_b)[14] = r2 >> 24;
           (p_b)[15] = r3 >> 24;
           (p_b)[16] = r0 >> 32;
           (p_b)[17] = r1 >> 32;
           (p_b)[18] = r2 >> 32;
           (p_b)[19] = r3 >> 32;
           (p_b)[20] = r0 >> 40;
           (p_b)[21] = r1 >> 40;
           (p_b)[22] = r2 >> 40;
           (p_b)[23] = r3 >> 40;
           p_src += 1;
           p_b += 32;
        }

        for (j = 0; j < 16; ++j) {
           // %%OP_XI_STORE b, p_perm, p_sign, p_dest, {i_dest.SHAPE[2]}
           uint_mmv_t r0, r1;
           r0 = ((uint_mmv_t)(b[p_perm[0]]) << 0)
             + ((uint_mmv_t)(b[p_perm[1]]) << 2)
             + ((uint_mmv_t)(b[p_perm[2]]) << 4)
             + ((uint_mmv_t)(b[p_perm[3]]) << 6)
             + ((uint_mmv_t)(b[p_perm[4]]) << 8)
             + ((uint_mmv_t)(b[p_perm[5]]) << 10)
             + ((uint_mmv_t)(b[p_perm[6]]) << 12)
             + ((uint_mmv_t)(b[p_perm[7]]) << 14)
             + ((uint_mmv_t)(b[p_perm[8]]) << 16)
             + ((uint_mmv_t)(b[p_perm[9]]) << 18)
             + ((uint_mmv_t)(b[p_perm[10]]) << 20)
             + ((uint_mmv_t)(b[p_perm[11]]) << 22)
             + ((uint_mmv_t)(b[p_perm[12]]) << 24)
             + ((uint_mmv_t)(b[p_perm[13]]) << 26)
             + ((uint_mmv_t)(b[p_perm[14]]) << 28)
             + ((uint_mmv_t)(b[p_perm[15]]) << 30)
             + ((uint_mmv_t)(b[p_perm[16]]) << 32)
             + ((uint_mmv_t)(b[p_perm[17]]) << 34)
             + ((uint_mmv_t)(b[p_perm[18]]) << 36)
             + ((uint_mmv_t)(b[p_perm[19]]) << 38)
             + ((uint_mmv_t)(b[p_perm[20]]) << 40)
             + ((uint_mmv_t)(b[p_perm[21]]) << 42)
             + ((uint_mmv_t)(b[p_perm[22]]) << 44)
             + ((uint_mmv_t)(b[p_perm[23]]) << 46);
           r1 = (p_sign)[0] >> 0;
           // Spread bits 0,...,31 of r1 to the (2-bit long) fields
           // of r1. A field of r1 is set to 0x3 if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xffffULL) 
               +  ((r1 & 0xffff0000ULL) << 16);
           r1 = (r1 & 0xff000000ffULL) 
               +  ((r1 & 0xff000000ff00ULL) << 8);
           r1 = (r1 & 0xf000f000f000fULL) 
               +  ((r1 & 0xf000f000f000f0ULL) << 4);
           r1 = (r1 & 0x303030303030303ULL) 
               +  ((r1 & 0xc0c0c0c0c0c0c0cULL) << 2);
           r1 = (r1 & 0x1111111111111111ULL) 
               +  ((r1 & 0x2222222222222222ULL) << 1);
           r1 = (((r1) << 2) - (r1));
           // Bit spreading done.
           (p_dest)[0] = r0 ^ r1;
           p_dest += 1;
           p_perm += 24;
           p_sign += 1;
        }
    }

    ///////////////////////////////////////////////////////////////
    // Map tag X1 to tag T1 if e = 1
    // Map tag X0 to tag T1 if e = 2
    ///////////////////////////////////////////////////////////////
    p_src = v_in + 2614;
    p_src -= diff;
    p_dest = v_out + 822;
    p_sign = p_tables[4].p_sign;
    p_perm = p_tables[4].p_perm;

    for (i = 0; i < 64; ++i) {
        p_b = b;
        for (j = 0; j < 16; ++j) {
           // %%OP_XI_LOAD p_src, p_b, {i_src.SHAPE[2]}
           uint_mmv_t r0, r1, r2, r3;
           r0 =  (p_src)[0];
           r1 = (r0 >> 2) & 0x303030303030303ULL;
           r2 = (r0 >> 4) & 0x303030303030303ULL;
           r3 = (r0 >> 6) & 0x303030303030303ULL;
           r0 &= 0x303030303030303ULL;
           (p_b)[0] = r0 >> 0;
           (p_b)[1] = r1 >> 0;
           (p_b)[2] = r2 >> 0;
           (p_b)[3] = r3 >> 0;
           (p_b)[4] = r0 >> 8;
           (p_b)[5] = r1 >> 8;
           (p_b)[6] = r2 >> 8;
           (p_b)[7] = r3 >> 8;
           (p_b)[8] = r0 >> 16;
           (p_b)[9] = r1 >> 16;
           (p_b)[10] = r2 >> 16;
           (p_b)[11] = r3 >> 16;
           (p_b)[12] = r0 >> 24;
           (p_b)[13] = r1 >> 24;
           (p_b)[14] = r2 >> 24;
           (p_b)[15] = r3 >> 24;
           (p_b)[16] = r0 >> 32;
           (p_b)[17] = r1 >> 32;
           (p_b)[18] = r2 >> 32;
           (p_b)[19] = r3 >> 32;
           (p_b)[20] = r0 >> 40;
           (p_b)[21] = r1 >> 40;
           (p_b)[22] = r2 >> 40;
           (p_b)[23] = r3 >> 40;
           p_src += 1;
           p_b += 32;
        }

        for (j = 0; j < 12; ++j) {
           // %%OP_XI_STORE b, p_perm, p_sign, p_dest, {i_dest.SHAPE[2]}
           uint_mmv_t r0, r1;
           r0 = ((uint_mmv_t)(b[p_perm[0]]) << 0)
             + ((uint_mmv_t)(b[p_perm[1]]) << 2)
             + ((uint_mmv_t)(b[p_perm[2]]) << 4)
             + ((uint_mmv_t)(b[p_perm[3]]) << 6)
             + ((uint_mmv_t)(b[p_perm[4]]) << 8)
             + ((uint_mmv_t)(b[p_perm[5]]) << 10)
             + ((uint_mmv_t)(b[p_perm[6]]) << 12)
             + ((uint_mmv_t)(b[p_perm[7]]) << 14)
             + ((uint_mmv_t)(b[p_perm[8]]) << 16)
             + ((uint_mmv_t)(b[p_perm[9]]) << 18)
             + ((uint_mmv_t)(b[p_perm[10]]) << 20)
             + ((uint_mmv_t)(b[p_perm[11]]) << 22)
             + ((uint_mmv_t)(b[p_perm[12]]) << 24)
             + ((uint_mmv_t)(b[p_perm[13]]) << 26)
             + ((uint_mmv_t)(b[p_perm[14]]) << 28)
             + ((uint_mmv_t)(b[p_perm[15]]) << 30)
             + ((uint_mmv_t)(b[p_perm[16]]) << 32)
             + ((uint_mmv_t)(b[p_perm[17]]) << 34)
             + ((uint_mmv_t)(b[p_perm[18]]) << 36)
             + ((uint_mmv_t)(b[p_perm[19]]) << 38)
             + ((uint_mmv_t)(b[p_perm[20]]) << 40)
             + ((uint_mmv_t)(b[p_perm[21]]) << 42)
             + ((uint_mmv_t)(b[p_perm[22]]) << 44)
             + ((uint_mmv_t)(b[p_perm[23]]) << 46)
             + ((uint_mmv_t)(b[p_perm[24]]) << 48)
             + ((uint_mmv_t)(b[p_perm[25]]) << 50)
             + ((uint_mmv_t)(b[p_perm[26]]) << 52)
             + ((uint_mmv_t)(b[p_perm[27]]) << 54)
             + ((uint_mmv_t)(b[p_perm[28]]) << 56)
             + ((uint_mmv_t)(b[p_perm[29]]) << 58)
             + ((uint_mmv_t)(b[p_perm[30]]) << 60)
             + ((uint_mmv_t)(b[p_perm[31]]) << 62);
           r1 = (p_sign)[0] >> 0;
           // Spread bits 0,...,31 of r1 to the (2-bit long) fields
           // of r1. A field of r1 is set to 0x3 if its 
           // corresponding bit in input r1 is one and to 0 otherwise.
           r1 = (r1 & 0xffffULL) 
               +  ((r1 & 0xffff0000ULL) << 16);
           r1 = (r1 & 0xff000000ffULL) 
               +  ((r1 & 0xff000000ff00ULL) << 8);
           r1 = (r1 & 0xf000f000f000fULL) 
               +  ((r1 & 0xf000f000f000f0ULL) << 4);
           r1 = (r1 & 0x303030303030303ULL) 
               +  ((r1 & 0xc0c0c0c0c0c0c0cULL) << 2);
           r1 = (r1 & 0x1111111111111111ULL) 
               +  ((r1 & 0x2222222222222222ULL) << 1);
           r1 = (((r1) << 2) - (r1));
           // Bit spreading done.
           (p_dest)[0] = r0 ^ r1;
           p_dest += 1;
           p_perm += 32;
           p_sign += 1;
        }
    }
}

static uint_mmv_t TAB3_XI64_MASK[] = {
// %%TABLE TABLE_MUL_MATRIX_XI64, uint{INT_BITS}
0x0000030303030303ULL,0x0000000000000000ULL,
0x0000030303030303ULL,0x0000ffffffffffffULL,
0x0000000000000000ULL,0x0000030303030303ULL,
0x0000ffffffffffffULL,0x0000030303030303ULL
};


#define HALF_YZ_SHIFT 10

static uint32_t TAB3_XI64_OFFSET[2][4] = {
    {
        MM_OP3_OFS_Z + (2 << HALF_YZ_SHIFT),
        MM_OP3_OFS_Z + (0 << HALF_YZ_SHIFT),
        MM_OP3_OFS_Z + (1 << HALF_YZ_SHIFT),
        MM_OP3_OFS_Z + (3 << HALF_YZ_SHIFT)
    },
    {
        MM_OP3_OFS_Z + (1 << HALF_YZ_SHIFT),
        MM_OP3_OFS_Z + (2 << HALF_YZ_SHIFT),
        MM_OP3_OFS_Z + (0 << HALF_YZ_SHIFT),
        MM_OP3_OFS_Z + (3 << HALF_YZ_SHIFT)
    },
};




static void mm_op3_xi_yz(uint_mmv_t *v_in,  uint32_t exp1, uint_mmv_t *v_out)
{
    uint_fast32_t i;
    uint_mmv_t *p_mask =  TAB3_XI64_MASK + exp1;
    for (i = 0; i < 64; ++i) {
        // %%MUL_MATRIX_XI64 v_in, p_mask, v_out

        // This is an automatically generated matrix operation, do not change!
        {
        uint_mmv_t r0, r1, r2, r3, r4;
        uint_mmv_t r5, r6, r7, r8, r9;

        uint_mmv_t a[24];
        uint_fast32_t i;

        // TODO: write comment!!!
        // 
        r0 = v_in[0] ^  p_mask[2];
        r9 = ((r0 ^ (r0 >> 2)) & 0xc0c0c0c0c0cULL);
        r0 ^= (r9 | (r9 << 2));
        r8 = v_in[14] ^  p_mask[0];
        r9 = ((r8 ^ (r8 >> 2)) & 0xc0c0c0c0c0cULL);
        r8 ^= (r9 | (r9 << 2));
        a[8] = (((r0 >> 16) & 0xffffULL)
            | ((r8 << 16) & 0xffff00000000ULL));
        a[16] = (((r0 >> 32) & 0xffffULL)
            | (r8 & 0xffff00000000ULL));
        r0 = ((r0 & 0xffffULL)
            | ((r8 << 32) & 0xffff00000000ULL));
        r1 = v_in[13] ^  p_mask[0];
        r9 = ((r1 ^ (r1 >> 2)) & 0xc0c0c0c0c0cULL);
        r1 ^= (r9 | (r9 << 2));
        r8 = v_in[3] ^  p_mask[0];
        r9 = ((r8 ^ (r8 >> 2)) & 0xc0c0c0c0c0cULL);
        r8 ^= (r9 | (r9 << 2));
        a[9] = (((r1 >> 16) & 0xffffULL)
            | ((r8 << 16) & 0xffff00000000ULL));
        a[17] = (((r1 >> 32) & 0xffffULL)
            | (r8 & 0xffff00000000ULL));
        r1 = ((r1 & 0xffffULL)
            | ((r8 << 32) & 0xffff00000000ULL));
        r2 = v_in[11] ^  p_mask[0];
        r9 = ((r2 ^ (r2 >> 2)) & 0xc0c0c0c0c0cULL);
        r2 ^= (r9 | (r9 << 2));
        r8 = v_in[5] ^  p_mask[0];
        r9 = ((r8 ^ (r8 >> 2)) & 0xc0c0c0c0c0cULL);
        r8 ^= (r9 | (r9 << 2));
        a[10] = (((r2 >> 16) & 0xffffULL)
            | ((r8 << 16) & 0xffff00000000ULL));
        a[18] = (((r2 >> 32) & 0xffffULL)
            | (r8 & 0xffff00000000ULL));
        r2 = ((r2 & 0xffffULL)
            | ((r8 << 32) & 0xffff00000000ULL));
        r3 = v_in[6] ^  p_mask[0];
        r9 = ((r3 ^ (r3 >> 2)) & 0xc0c0c0c0c0cULL);
        r3 ^= (r9 | (r9 << 2));
        r8 = v_in[8] ^  p_mask[2];
        r9 = ((r8 ^ (r8 >> 2)) & 0xc0c0c0c0c0cULL);
        r8 ^= (r9 | (r9 << 2));
        a[11] = (((r3 >> 16) & 0xffffULL)
            | ((r8 << 16) & 0xffff00000000ULL));
        a[19] = (((r3 >> 32) & 0xffffULL)
            | (r8 & 0xffff00000000ULL));
        r3 = ((r3 & 0xffffULL)
            | ((r8 << 32) & 0xffff00000000ULL));
        r4 = v_in[7] ^  p_mask[0];
        r9 = ((r4 ^ (r4 >> 2)) & 0xc0c0c0c0c0cULL);
        r4 ^= (r9 | (r9 << 2));
        r8 = v_in[9] ^  p_mask[0];
        r9 = ((r8 ^ (r8 >> 2)) & 0xc0c0c0c0c0cULL);
        r8 ^= (r9 | (r9 << 2));
        a[12] = (((r4 >> 16) & 0xffffULL)
            | ((r8 << 16) & 0xffff00000000ULL));
        a[20] = (((r4 >> 32) & 0xffffULL)
            | (r8 & 0xffff00000000ULL));
        r4 = ((r4 & 0xffffULL)
            | ((r8 << 32) & 0xffff00000000ULL));
        r5 = v_in[10] ^  p_mask[0];
        r9 = ((r5 ^ (r5 >> 2)) & 0xc0c0c0c0c0cULL);
        r5 ^= (r9 | (r9 << 2));
        r8 = v_in[4] ^  p_mask[2];
        r9 = ((r8 ^ (r8 >> 2)) & 0xc0c0c0c0c0cULL);
        r8 ^= (r9 | (r9 << 2));
        a[13] = (((r5 >> 16) & 0xffffULL)
            | ((r8 << 16) & 0xffff00000000ULL));
        a[21] = (((r5 >> 32) & 0xffffULL)
            | (r8 & 0xffff00000000ULL));
        r5 = ((r5 & 0xffffULL)
            | ((r8 << 32) & 0xffff00000000ULL));
        r6 = v_in[12] ^  p_mask[0];
        r9 = ((r6 ^ (r6 >> 2)) & 0xc0c0c0c0c0cULL);
        r6 ^= (r9 | (r9 << 2));
        r8 = v_in[2] ^  p_mask[2];
        r9 = ((r8 ^ (r8 >> 2)) & 0xc0c0c0c0c0cULL);
        r8 ^= (r9 | (r9 << 2));
        a[14] = (((r6 >> 16) & 0xffffULL)
            | ((r8 << 16) & 0xffff00000000ULL));
        a[22] = (((r6 >> 32) & 0xffffULL)
            | (r8 & 0xffff00000000ULL));
        r6 = ((r6 & 0xffffULL)
            | ((r8 << 32) & 0xffff00000000ULL));
        r7 = v_in[1] ^  p_mask[2];
        r9 = ((r7 ^ (r7 >> 2)) & 0xc0c0c0c0c0cULL);
        r7 ^= (r9 | (r9 << 2));
        r8 = v_in[15] ^  p_mask[2];
        r9 = ((r8 ^ (r8 >> 2)) & 0xc0c0c0c0c0cULL);
        r8 ^= (r9 | (r9 << 2));
        a[15] = (((r7 >> 16) & 0xffffULL)
            | ((r8 << 16) & 0xffff00000000ULL));
        a[23] = (((r7 >> 32) & 0xffffULL)
            | (r8 & 0xffff00000000ULL));
        r7 = ((r7 & 0xffffULL)
            | ((r8 << 32) & 0xffff00000000ULL));
        // 72 lines of code, 216 operations
        i = 0;
        goto l_mmv3_op_l64_2;
        l_mmv3_op_l64_1:
        a[(i) + 0] = r0;
        a[(i) + 1] = r1;
        a[(i) + 2] = r2;
        a[(i) + 3] = r3;
        a[(i) + 4] = r4;
        a[(i) + 5] = r5;
        a[(i) + 6] = r6;
        a[(i) + 7] = r7;
        i += 8;
        r0 = a[(i) + 0];
        r1 = a[(i) + 1];
        r2 = a[(i) + 2];
        r3 = a[(i) + 3];
        r4 = a[(i) + 4];
        r5 = a[(i) + 5];
        r6 = a[(i) + 6];
        r7 = a[(i) + 7];
        l_mmv3_op_l64_2:
        // Expansion for Hadamard operation:
        // There is no space for a carry bit between bit fields. So 
        // we move bit field 2*i + 1  to bit field 2*i + 8.
        r0 = ((r0 & 0x333300003333ULL)
            | ((r0 & 0xcccc0000ccccULL) << 14));
        r1 = ((r1 & 0x333300003333ULL)
            | ((r1 & 0xcccc0000ccccULL) << 14));
        r2 = ((r2 & 0x333300003333ULL)
            | ((r2 & 0xcccc0000ccccULL) << 14));
        r3 = ((r3 & 0x333300003333ULL)
            | ((r3 & 0xcccc0000ccccULL) << 14));
        r4 = ((r4 & 0x333300003333ULL)
            | ((r4 & 0xcccc0000ccccULL) << 14));
        r5 = ((r5 & 0x333300003333ULL)
            | ((r5 & 0xcccc0000ccccULL) << 14));
        r6 = ((r6 & 0x333300003333ULL)
            | ((r6 & 0xcccc0000ccccULL) << 14));
        r7 = ((r7 & 0x333300003333ULL)
            | ((r7 & 0xcccc0000ccccULL) << 14));
        // Vector is now  r(i) for i = 0,1,2,3,4,5,6,7
        // Butterfly: v[i], v[i+2] = v[i]+v[i+2], v[i]-v[i+2]
        r8 = (((r0 << 4) & 0x3030303030303030ULL)
            | ((r0 & 0x3030303030303030ULL) >> 4));
        r0 = ((r0 ^ 0x3030303030303030ULL) + r8);
        r8 = (r0 & 0x4444444444444444ULL);
        r0 = ((r0 - r8) + (r8 >> 2));
        r8 = (((r1 << 4) & 0x3030303030303030ULL)
            | ((r1 & 0x3030303030303030ULL) >> 4));
        r1 = ((r1 ^ 0x3030303030303030ULL) + r8);
        r8 = (r1 & 0x4444444444444444ULL);
        r1 = ((r1 - r8) + (r8 >> 2));
        r8 = (((r2 << 4) & 0x3030303030303030ULL)
            | ((r2 & 0x3030303030303030ULL) >> 4));
        r2 = ((r2 ^ 0x3030303030303030ULL) + r8);
        r8 = (r2 & 0x4444444444444444ULL);
        r2 = ((r2 - r8) + (r8 >> 2));
        r8 = (((r3 << 4) & 0x3030303030303030ULL)
            | ((r3 & 0x3030303030303030ULL) >> 4));
        r3 = ((r3 ^ 0x3030303030303030ULL) + r8);
        r8 = (r3 & 0x4444444444444444ULL);
        r3 = ((r3 - r8) + (r8 >> 2));
        r8 = (((r4 << 4) & 0x3030303030303030ULL)
            | ((r4 & 0x3030303030303030ULL) >> 4));
        r4 = ((r4 ^ 0x3030303030303030ULL) + r8);
        r8 = (r4 & 0x4444444444444444ULL);
        r4 = ((r4 - r8) + (r8 >> 2));
        r8 = (((r5 << 4) & 0x3030303030303030ULL)
            | ((r5 & 0x3030303030303030ULL) >> 4));
        r5 = ((r5 ^ 0x3030303030303030ULL) + r8);
        r8 = (r5 & 0x4444444444444444ULL);
        r5 = ((r5 - r8) + (r8 >> 2));
        r8 = (((r6 << 4) & 0x3030303030303030ULL)
            | ((r6 & 0x3030303030303030ULL) >> 4));
        r6 = ((r6 ^ 0x3030303030303030ULL) + r8);
        r8 = (r6 & 0x4444444444444444ULL);
        r6 = ((r6 - r8) + (r8 >> 2));
        r8 = (((r7 << 4) & 0x3030303030303030ULL)
            | ((r7 & 0x3030303030303030ULL) >> 4));
        r7 = ((r7 ^ 0x3030303030303030ULL) + r8);
        r8 = (r7 & 0x4444444444444444ULL);
        r7 = ((r7 - r8) + (r8 >> 2));
        // Vector is now  r(i) for i = 0,1,2,3,4,5,6,7
        // Butterfly: v[i], v[i+8] = v[i]+v[i+8], v[i]-v[i+8]
        r8 = (((r0 << 16) & 0x3333000033330000ULL)
            | ((r0 & 0x3333000033330000ULL) >> 16));
        r0 = ((r0 ^ 0x3333000033330000ULL) + r8);
        r8 = (r0 & 0x4444444444444444ULL);
        r0 = ((r0 - r8) + (r8 >> 2));
        r8 = (((r1 << 16) & 0x3333000033330000ULL)
            | ((r1 & 0x3333000033330000ULL) >> 16));
        r1 = ((r1 ^ 0x3333000033330000ULL) + r8);
        r8 = (r1 & 0x4444444444444444ULL);
        r1 = ((r1 - r8) + (r8 >> 2));
        r8 = (((r2 << 16) & 0x3333000033330000ULL)
            | ((r2 & 0x3333000033330000ULL) >> 16));
        r2 = ((r2 ^ 0x3333000033330000ULL) + r8);
        r8 = (r2 & 0x4444444444444444ULL);
        r2 = ((r2 - r8) + (r8 >> 2));
        r8 = (((r3 << 16) & 0x3333000033330000ULL)
            | ((r3 & 0x3333000033330000ULL) >> 16));
        r3 = ((r3 ^ 0x3333000033330000ULL) + r8);
        r8 = (r3 & 0x4444444444444444ULL);
        r3 = ((r3 - r8) + (r8 >> 2));
        r8 = (((r4 << 16) & 0x3333000033330000ULL)
            | ((r4 & 0x3333000033330000ULL) >> 16));
        r4 = ((r4 ^ 0x3333000033330000ULL) + r8);
        r8 = (r4 & 0x4444444444444444ULL);
        r4 = ((r4 - r8) + (r8 >> 2));
        r8 = (((r5 << 16) & 0x3333000033330000ULL)
            | ((r5 & 0x3333000033330000ULL) >> 16));
        r5 = ((r5 ^ 0x3333000033330000ULL) + r8);
        r8 = (r5 & 0x4444444444444444ULL);
        r5 = ((r5 - r8) + (r8 >> 2));
        r8 = (((r6 << 16) & 0x3333000033330000ULL)
            | ((r6 & 0x3333000033330000ULL) >> 16));
        r6 = ((r6 ^ 0x3333000033330000ULL) + r8);
        r8 = (r6 & 0x4444444444444444ULL);
        r6 = ((r6 - r8) + (r8 >> 2));
        r8 = (((r7 << 16) & 0x3333000033330000ULL)
            | ((r7 & 0x3333000033330000ULL) >> 16));
        r7 = ((r7 ^ 0x3333000033330000ULL) + r8);
        r8 = (r7 & 0x4444444444444444ULL);
        r7 = ((r7 - r8) + (r8 >> 2));
        // Vector is now  r(i) for i = 0,1,2,3,4,5,6,7
        // Butterfly: v[i], v[i+16] = v[i]+v[i+16], v[i]-v[i+16]
        r8 = ((r0 << 32) | (r0 >> 32));
        r0 = ((r0 ^ 0x3333333300000000ULL) + r8);
        r8 = (r0 & 0x4444444444444444ULL);
        r0 = ((r0 - r8) + (r8 >> 2));
        r8 = ((r1 << 32) | (r1 >> 32));
        r1 = ((r1 ^ 0x3333333300000000ULL) + r8);
        r8 = (r1 & 0x4444444444444444ULL);
        r1 = ((r1 - r8) + (r8 >> 2));
        r8 = ((r2 << 32) | (r2 >> 32));
        r2 = ((r2 ^ 0x3333333300000000ULL) + r8);
        r8 = (r2 & 0x4444444444444444ULL);
        r2 = ((r2 - r8) + (r8 >> 2));
        r8 = ((r3 << 32) | (r3 >> 32));
        r3 = ((r3 ^ 0x3333333300000000ULL) + r8);
        r8 = (r3 & 0x4444444444444444ULL);
        r3 = ((r3 - r8) + (r8 >> 2));
        r8 = ((r4 << 32) | (r4 >> 32));
        r4 = ((r4 ^ 0x3333333300000000ULL) + r8);
        r8 = (r4 & 0x4444444444444444ULL);
        r4 = ((r4 - r8) + (r8 >> 2));
        r8 = ((r5 << 32) | (r5 >> 32));
        r5 = ((r5 ^ 0x3333333300000000ULL) + r8);
        r8 = (r5 & 0x4444444444444444ULL);
        r5 = ((r5 - r8) + (r8 >> 2));
        r8 = ((r6 << 32) | (r6 >> 32));
        r6 = ((r6 ^ 0x3333333300000000ULL) + r8);
        r8 = (r6 & 0x4444444444444444ULL);
        r6 = ((r6 - r8) + (r8 >> 2));
        r8 = ((r7 << 32) | (r7 >> 32));
        r7 = ((r7 ^ 0x3333333300000000ULL) + r8);
        r8 = (r7 & 0x4444444444444444ULL);
        r7 = ((r7 - r8) + (r8 >> 2));
        // Vector is now  r(i) for i = 0,1,2,3,4,5,6,7
        // Butterfly: v[i], v[i+32] = v[i]+v[i+32], v[i]-v[i+32]
        r8 = (r0 + (r1 ^ 0x3333333333333333ULL));
        r0 = (r0 + r1);
        r1 = (r8 & 0x4444444444444444ULL);
        r1 = ((r8 - r1) + (r1 >> 2));
        r8 = (r0 & 0x4444444444444444ULL);
        r0 = ((r0 - r8) + (r8 >> 2));
        r8 = (r2 + (r3 ^ 0x3333333333333333ULL));
        r2 = (r2 + r3);
        r3 = (r8 & 0x4444444444444444ULL);
        r3 = ((r8 - r3) + (r3 >> 2));
        r8 = (r2 & 0x4444444444444444ULL);
        r2 = ((r2 - r8) + (r8 >> 2));
        r8 = (r4 + (r5 ^ 0x3333333333333333ULL));
        r4 = (r4 + r5);
        r5 = (r8 & 0x4444444444444444ULL);
        r5 = ((r8 - r5) + (r5 >> 2));
        r8 = (r4 & 0x4444444444444444ULL);
        r4 = ((r4 - r8) + (r8 >> 2));
        r8 = (r6 + (r7 ^ 0x3333333333333333ULL));
        r6 = (r6 + r7);
        r7 = (r8 & 0x4444444444444444ULL);
        r7 = ((r8 - r7) + (r7 >> 2));
        r8 = (r6 & 0x4444444444444444ULL);
        r6 = ((r6 - r8) + (r8 >> 2));
        // Vector is now  r(i) for i = 0,1,2,3,4,5,6,7
        // Butterfly: v[i], v[i+64] = v[i]+v[i+64], v[i]-v[i+64]
        r8 = (r0 + (r2 ^ 0x3333333333333333ULL));
        r0 = (r0 + r2);
        r2 = (r8 & 0x4444444444444444ULL);
        r2 = ((r8 - r2) + (r2 >> 2));
        r8 = (r0 & 0x4444444444444444ULL);
        r0 = ((r0 - r8) + (r8 >> 2));
        r8 = (r1 + (r3 ^ 0x3333333333333333ULL));
        r1 = (r1 + r3);
        r3 = (r8 & 0x4444444444444444ULL);
        r3 = ((r8 - r3) + (r3 >> 2));
        r8 = (r1 & 0x4444444444444444ULL);
        r1 = ((r1 - r8) + (r8 >> 2));
        r8 = (r4 + (r6 ^ 0x3333333333333333ULL));
        r4 = (r4 + r6);
        r6 = (r8 & 0x4444444444444444ULL);
        r6 = ((r8 - r6) + (r6 >> 2));
        r8 = (r4 & 0x4444444444444444ULL);
        r4 = ((r4 - r8) + (r8 >> 2));
        r8 = (r5 + (r7 ^ 0x3333333333333333ULL));
        r5 = (r5 + r7);
        r7 = (r8 & 0x4444444444444444ULL);
        r7 = ((r8 - r7) + (r7 >> 2));
        r8 = (r5 & 0x4444444444444444ULL);
        r5 = ((r5 - r8) + (r8 >> 2));
        // Vector is now  r(i) for i = 0,1,2,3,4,5,6,7
        // Butterfly: v[i], v[i+128] = v[i]+v[i+128], v[i]-v[i+128]
        r8 = (r0 + (r4 ^ 0x3333333333333333ULL));
        r0 = (r0 + r4);
        r4 = (r8 & 0x4444444444444444ULL);
        r4 = ((r8 - r4) + (r4 >> 2));
        r8 = (r0 & 0x4444444444444444ULL);
        r0 = ((r0 - r8) + (r8 >> 2));
        r8 = (r1 + (r5 ^ 0x3333333333333333ULL));
        r1 = (r1 + r5);
        r5 = (r8 & 0x4444444444444444ULL);
        r5 = ((r8 - r5) + (r5 >> 2));
        r8 = (r1 & 0x4444444444444444ULL);
        r1 = ((r1 - r8) + (r8 >> 2));
        r8 = (r2 + (r6 ^ 0x3333333333333333ULL));
        r2 = (r2 + r6);
        r6 = (r8 & 0x4444444444444444ULL);
        r6 = ((r8 - r6) + (r6 >> 2));
        r8 = (r2 & 0x4444444444444444ULL);
        r2 = ((r2 - r8) + (r8 >> 2));
        r8 = (r3 + (r7 ^ 0x3333333333333333ULL));
        r3 = (r3 + r7);
        r7 = (r8 & 0x4444444444444444ULL);
        r7 = ((r8 - r7) + (r7 >> 2));
        r8 = (r3 & 0x4444444444444444ULL);
        r3 = ((r3 - r8) + (r8 >> 2));
        // Vector is now  r(i) for i = 0,1,2,3,4,5,6,7
        // Reverse expansion for Hadamard operation
        r0 = ((r0 & 0x333300003333ULL)
            | ((r0 >> 14) & 0xcccc0000ccccULL));
        r1 = ((r1 & 0x333300003333ULL)
            | ((r1 >> 14) & 0xcccc0000ccccULL));
        r2 = ((r2 & 0x333300003333ULL)
            | ((r2 >> 14) & 0xcccc0000ccccULL));
        r3 = ((r3 & 0x333300003333ULL)
            | ((r3 >> 14) & 0xcccc0000ccccULL));
        r4 = ((r4 & 0x333300003333ULL)
            | ((r4 >> 14) & 0xcccc0000ccccULL));
        r5 = ((r5 & 0x333300003333ULL)
            | ((r5 >> 14) & 0xcccc0000ccccULL));
        r6 = ((r6 & 0x333300003333ULL)
            | ((r6 >> 14) & 0xcccc0000ccccULL));
        r7 = ((r7 & 0x333300003333ULL)
            | ((r7 >> 14) & 0xcccc0000ccccULL));
        // Vector is now  r(i) for i = 0,1,2,3,4,5,6,7
        // 184 lines of code, 444 operations
        if (i < 16) goto l_mmv3_op_l64_1;
        r8 = (((a[0] & 0xffffULL)
            | ((a[8] << 16) & 0xffff0000ULL))
            | (r0 << 32));
        r8 = (((r8 & 0x5555555555555555ULL) << 1)
            | ((r8 & 0xaaaaaaaaaaaaaaaaULL) >> 1));
        v_out[0] = (r8 ^  p_mask[6]) & 0xffffffffffffULL;
        a[0] = ((((a[0] >> 32) & 0xffffULL)
            | ((a[8] >> 16) & 0xffff0000ULL))
            | (r0 & 0xffff00000000ULL));
        a[0] = (((a[0] & 0x5555555555555555ULL) << 1)
            | ((a[0] & 0xaaaaaaaaaaaaaaaaULL) >> 1));
        v_out[1] = (a[0] ^  p_mask[6]) & 0xffffffffffffULL;
        r8 = (((a[1] & 0xffffULL)
            | ((a[9] << 16) & 0xffff0000ULL))
            | (r1 << 32));
        r8 = (((r8 & 0x5555555555555555ULL) << 1)
            | ((r8 & 0xaaaaaaaaaaaaaaaaULL) >> 1));
        v_out[2] = (r8 ^  p_mask[6]) & 0xffffffffffffULL;
        a[1] = ((((a[1] >> 32) & 0xffffULL)
            | ((a[9] >> 16) & 0xffff0000ULL))
            | (r1 & 0xffff00000000ULL));
        a[1] = (((a[1] & 0x5555555555555555ULL) << 1)
            | ((a[1] & 0xaaaaaaaaaaaaaaaaULL) >> 1));
        v_out[3] = (a[1] ^  p_mask[4]) & 0xffffffffffffULL;
        r8 = (((a[2] & 0xffffULL)
            | ((a[10] << 16) & 0xffff0000ULL))
            | (r2 << 32));
        r8 = (((r8 & 0x5555555555555555ULL) << 1)
            | ((r8 & 0xaaaaaaaaaaaaaaaaULL) >> 1));
        v_out[4] = (r8 ^  p_mask[6]) & 0xffffffffffffULL;
        a[2] = ((((a[2] >> 32) & 0xffffULL)
            | ((a[10] >> 16) & 0xffff0000ULL))
            | (r2 & 0xffff00000000ULL));
        a[2] = (((a[2] & 0x5555555555555555ULL) << 1)
            | ((a[2] & 0xaaaaaaaaaaaaaaaaULL) >> 1));
        v_out[5] = (a[2] ^  p_mask[4]) & 0xffffffffffffULL;
        r8 = (((a[3] & 0xffffULL)
            | ((a[11] << 16) & 0xffff0000ULL))
            | (r3 << 32));
        r8 = (((r8 & 0x5555555555555555ULL) << 1)
            | ((r8 & 0xaaaaaaaaaaaaaaaaULL) >> 1));
        v_out[6] = (r8 ^  p_mask[4]) & 0xffffffffffffULL;
        a[3] = ((((a[3] >> 32) & 0xffffULL)
            | ((a[11] >> 16) & 0xffff0000ULL))
            | (r3 & 0xffff00000000ULL));
        a[3] = (((a[3] & 0x5555555555555555ULL) << 1)
            | ((a[3] & 0xaaaaaaaaaaaaaaaaULL) >> 1));
        v_out[7] = (a[3] ^  p_mask[4]) & 0xffffffffffffULL;
        r8 = (((a[4] & 0xffffULL)
            | ((a[12] << 16) & 0xffff0000ULL))
            | (r4 << 32));
        r8 = (((r8 & 0x5555555555555555ULL) << 1)
            | ((r8 & 0xaaaaaaaaaaaaaaaaULL) >> 1));
        v_out[8] = (r8 ^  p_mask[6]) & 0xffffffffffffULL;
        a[4] = ((((a[4] >> 32) & 0xffffULL)
            | ((a[12] >> 16) & 0xffff0000ULL))
            | (r4 & 0xffff00000000ULL));
        a[4] = (((a[4] & 0x5555555555555555ULL) << 1)
            | ((a[4] & 0xaaaaaaaaaaaaaaaaULL) >> 1));
        v_out[9] = (a[4] ^  p_mask[4]) & 0xffffffffffffULL;
        r8 = (((a[5] & 0xffffULL)
            | ((a[13] << 16) & 0xffff0000ULL))
            | (r5 << 32));
        r8 = (((r8 & 0x5555555555555555ULL) << 1)
            | ((r8 & 0xaaaaaaaaaaaaaaaaULL) >> 1));
        v_out[10] = (r8 ^  p_mask[4]) & 0xffffffffffffULL;
        a[5] = ((((a[5] >> 32) & 0xffffULL)
            | ((a[13] >> 16) & 0xffff0000ULL))
            | (r5 & 0xffff00000000ULL));
        a[5] = (((a[5] & 0x5555555555555555ULL) << 1)
            | ((a[5] & 0xaaaaaaaaaaaaaaaaULL) >> 1));
        v_out[11] = (a[5] ^  p_mask[4]) & 0xffffffffffffULL;
        r8 = (((a[6] & 0xffffULL)
            | ((a[14] << 16) & 0xffff0000ULL))
            | (r6 << 32));
        r8 = (((r8 & 0x5555555555555555ULL) << 1)
            | ((r8 & 0xaaaaaaaaaaaaaaaaULL) >> 1));
        v_out[12] = (r8 ^  p_mask[4]) & 0xffffffffffffULL;
        a[6] = ((((a[6] >> 32) & 0xffffULL)
            | ((a[14] >> 16) & 0xffff0000ULL))
            | (r6 & 0xffff00000000ULL));
        a[6] = (((a[6] & 0x5555555555555555ULL) << 1)
            | ((a[6] & 0xaaaaaaaaaaaaaaaaULL) >> 1));
        v_out[13] = (a[6] ^  p_mask[4]) & 0xffffffffffffULL;
        r8 = (((a[7] & 0xffffULL)
            | ((a[15] << 16) & 0xffff0000ULL))
            | (r7 << 32));
        r8 = (((r8 & 0x5555555555555555ULL) << 1)
            | ((r8 & 0xaaaaaaaaaaaaaaaaULL) >> 1));
        v_out[14] = (r8 ^  p_mask[4]) & 0xffffffffffffULL;
        a[7] = ((((a[7] >> 32) & 0xffffULL)
            | ((a[15] >> 16) & 0xffff0000ULL))
            | (r7 & 0xffff00000000ULL));
        a[7] = (((a[7] & 0x5555555555555555ULL) << 1)
            | ((a[7] & 0xaaaaaaaaaaaaaaaaULL) >> 1));
        v_out[15] = (a[7] ^  p_mask[6]) & 0xffffffffffffULL;
        v_in += 16;
        v_out += 16;
        // 48 lines of code, 224 operations
        }
        // End of automatically generated matrix operation.
 
    }
}  


static void mm_op3_xi_a(uint_mmv_t *v_in,  uint32_t exp1, uint_mmv_t *v_out)
{
    uint_fast32_t i;
    uint_mmv_t e_mask =  -((uint_mmv_t)exp1 & 0x1ULL);
    for (i = 0; i < 6; ++i) {
        // %%MUL_MATRIX_XI16 v_in, e_mask, v_out

        // This is an automatically generated matrix operation, do not change!
        {
        uint_mmv_t r0, r1, r2, r3, r4;
        uint_mmv_t r5, r6, r7, r8;

        // TODO: write comment!!!
        // 
        e_mask = ~(e_mask);
        r0 = v_in[0] ^  (0xfcfcfcfcfcfcULL & e_mask);
        r4 = ((r0 ^ (r0 >> 2)) & 0xc0c0c0c0c0cULL);
        r0 ^= (r4 | (r4 << 2));
        r1 = v_in[2] ^  (0x30303030303ULL & e_mask);
        r4 = ((r1 ^ (r1 >> 2)) & 0xc0c0c0c0c0cULL);
        r1 ^= (r4 | (r4 << 2));
        r2 = v_in[1] ^  (0x30303030303ULL & e_mask);
        r4 = ((r2 ^ (r2 >> 2)) & 0xc0c0c0c0c0cULL);
        r2 ^= (r4 | (r4 << 2));
        r3 = v_in[3] ^  (0x30303030303ULL & e_mask);
        r4 = ((r3 ^ (r3 >> 2)) & 0xc0c0c0c0c0cULL);
        r3 ^= (r4 | (r4 << 2));
        // Expansion for Hadamard operation:
        // There is no space for a carry bit between bit fields. So 
        // we move bit field 2*i + 1  to bit field 2*i + 128.
        r4 = ((r0 >> 2) & 0x3333333333333333ULL);
        r0 = (r0 & 0x3333333333333333ULL);
        r5 = ((r1 >> 2) & 0x3333333333333333ULL);
        r1 = (r1 & 0x3333333333333333ULL);
        r6 = ((r2 >> 2) & 0x3333333333333333ULL);
        r2 = (r2 & 0x3333333333333333ULL);
        r7 = ((r3 >> 2) & 0x3333333333333333ULL);
        r3 = (r3 & 0x3333333333333333ULL);
        // Vector is now  r(i) for i = 0,1,2,3,4,5,6,7
        // Butterfly: v[i], v[i+2] = v[i]+v[i+2], v[i]-v[i+2]
        r8 = (((r0 << 4) & 0x3030303030303030ULL)
            | ((r0 & 0x3030303030303030ULL) >> 4));
        r0 = ((r0 ^ 0x3030303030303030ULL) + r8);
        r8 = (r0 & 0x4444444444444444ULL);
        r0 = ((r0 - r8) + (r8 >> 2));
        r8 = (((r1 << 4) & 0x3030303030303030ULL)
            | ((r1 & 0x3030303030303030ULL) >> 4));
        r1 = ((r1 ^ 0x3030303030303030ULL) + r8);
        r8 = (r1 & 0x4444444444444444ULL);
        r1 = ((r1 - r8) + (r8 >> 2));
        r8 = (((r2 << 4) & 0x3030303030303030ULL)
            | ((r2 & 0x3030303030303030ULL) >> 4));
        r2 = ((r2 ^ 0x3030303030303030ULL) + r8);
        r8 = (r2 & 0x4444444444444444ULL);
        r2 = ((r2 - r8) + (r8 >> 2));
        r8 = (((r3 << 4) & 0x3030303030303030ULL)
            | ((r3 & 0x3030303030303030ULL) >> 4));
        r3 = ((r3 ^ 0x3030303030303030ULL) + r8);
        r8 = (r3 & 0x4444444444444444ULL);
        r3 = ((r3 - r8) + (r8 >> 2));
        r8 = (((r4 << 4) & 0x3030303030303030ULL)
            | ((r4 & 0x3030303030303030ULL) >> 4));
        r4 = ((r4 ^ 0x3030303030303030ULL) + r8);
        r8 = (r4 & 0x4444444444444444ULL);
        r4 = ((r4 - r8) + (r8 >> 2));
        r8 = (((r5 << 4) & 0x3030303030303030ULL)
            | ((r5 & 0x3030303030303030ULL) >> 4));
        r5 = ((r5 ^ 0x3030303030303030ULL) + r8);
        r8 = (r5 & 0x4444444444444444ULL);
        r5 = ((r5 - r8) + (r8 >> 2));
        r8 = (((r6 << 4) & 0x3030303030303030ULL)
            | ((r6 & 0x3030303030303030ULL) >> 4));
        r6 = ((r6 ^ 0x3030303030303030ULL) + r8);
        r8 = (r6 & 0x4444444444444444ULL);
        r6 = ((r6 - r8) + (r8 >> 2));
        r8 = (((r7 << 4) & 0x3030303030303030ULL)
            | ((r7 & 0x3030303030303030ULL) >> 4));
        r7 = ((r7 ^ 0x3030303030303030ULL) + r8);
        r8 = (r7 & 0x4444444444444444ULL);
        r7 = ((r7 - r8) + (r8 >> 2));
        // Vector is now  r(i) for i = 0,1,2,3,4,5,6,7
        // Butterfly: v[i], v[i+32] = v[i]+v[i+32], v[i]-v[i+32]
        r8 = (r0 + (r1 ^ 0x3333333333333333ULL));
        r0 = (r0 + r1);
        r1 = (r8 & 0x4444444444444444ULL);
        r1 = ((r8 - r1) + (r1 >> 2));
        r8 = (r0 & 0x4444444444444444ULL);
        r0 = ((r0 - r8) + (r8 >> 2));
        r8 = (r2 + (r3 ^ 0x3333333333333333ULL));
        r2 = (r2 + r3);
        r3 = (r8 & 0x4444444444444444ULL);
        r3 = ((r8 - r3) + (r3 >> 2));
        r8 = (r2 & 0x4444444444444444ULL);
        r2 = ((r2 - r8) + (r8 >> 2));
        r8 = (r4 + (r5 ^ 0x3333333333333333ULL));
        r4 = (r4 + r5);
        r5 = (r8 & 0x4444444444444444ULL);
        r5 = ((r8 - r5) + (r5 >> 2));
        r8 = (r4 & 0x4444444444444444ULL);
        r4 = ((r4 - r8) + (r8 >> 2));
        r8 = (r6 + (r7 ^ 0x3333333333333333ULL));
        r6 = (r6 + r7);
        r7 = (r8 & 0x4444444444444444ULL);
        r7 = ((r8 - r7) + (r7 >> 2));
        r8 = (r6 & 0x4444444444444444ULL);
        r6 = ((r6 - r8) + (r8 >> 2));
        // Vector is now  r(i) for i = 0,1,2,3,4,5,6,7
        // Butterfly: v[i], v[i+64] = v[i]+v[i+64], v[i]-v[i+64]
        r8 = (r0 + (r2 ^ 0x3333333333333333ULL));
        r0 = (r0 + r2);
        r2 = (r8 & 0x4444444444444444ULL);
        r2 = ((r8 - r2) + (r2 >> 2));
        r8 = (r0 & 0x4444444444444444ULL);
        r0 = ((r0 - r8) + (r8 >> 2));
        r8 = (r1 + (r3 ^ 0x3333333333333333ULL));
        r1 = (r1 + r3);
        r3 = (r8 & 0x4444444444444444ULL);
        r3 = ((r8 - r3) + (r3 >> 2));
        r8 = (r1 & 0x4444444444444444ULL);
        r1 = ((r1 - r8) + (r8 >> 2));
        r8 = (r4 + (r6 ^ 0x3333333333333333ULL));
        r4 = (r4 + r6);
        r6 = (r8 & 0x4444444444444444ULL);
        r6 = ((r8 - r6) + (r6 >> 2));
        r8 = (r4 & 0x4444444444444444ULL);
        r4 = ((r4 - r8) + (r8 >> 2));
        r8 = (r5 + (r7 ^ 0x3333333333333333ULL));
        r5 = (r5 + r7);
        r7 = (r8 & 0x4444444444444444ULL);
        r7 = ((r8 - r7) + (r7 >> 2));
        r8 = (r5 & 0x4444444444444444ULL);
        r5 = ((r5 - r8) + (r8 >> 2));
        // Vector is now  r(i) for i = 0,1,2,3,4,5,6,7
        // Butterfly: v[i], v[i+128] = v[i]+v[i+128], v[i]-v[i+128]
        r8 = (r0 + (r4 ^ 0x3333333333333333ULL));
        r0 = (r0 + r4);
        r4 = (r8 & 0x4444444444444444ULL);
        r4 = ((r8 - r4) + (r4 >> 2));
        r8 = (r0 & 0x4444444444444444ULL);
        r0 = ((r0 - r8) + (r8 >> 2));
        r8 = (r1 + (r5 ^ 0x3333333333333333ULL));
        r1 = (r1 + r5);
        r5 = (r8 & 0x4444444444444444ULL);
        r5 = ((r8 - r5) + (r5 >> 2));
        r8 = (r1 & 0x4444444444444444ULL);
        r1 = ((r1 - r8) + (r8 >> 2));
        r8 = (r2 + (r6 ^ 0x3333333333333333ULL));
        r2 = (r2 + r6);
        r6 = (r8 & 0x4444444444444444ULL);
        r6 = ((r8 - r6) + (r6 >> 2));
        r8 = (r2 & 0x4444444444444444ULL);
        r2 = ((r2 - r8) + (r8 >> 2));
        r8 = (r3 + (r7 ^ 0x3333333333333333ULL));
        r3 = (r3 + r7);
        r7 = (r8 & 0x4444444444444444ULL);
        r7 = ((r8 - r7) + (r7 >> 2));
        r8 = (r3 & 0x4444444444444444ULL);
        r3 = ((r3 - r8) + (r8 >> 2));
        // Vector is now  r(i) for i = 0,1,2,3,4,5,6,7
        // Reverse expansion for Hadamard operation
        r0 ^= (r4 << 2);
        r1 ^= (r5 << 2);
        r2 ^= (r6 << 2);
        r3 ^= (r7 << 2);
        // Vector is now  r(i) for i = 0,1,2,3
        e_mask = ~(e_mask);
        v_out[0] = r0 ^ (e_mask & 0xfcfcfcfcfcfcULL);
        v_out[1] = r1 ^ (e_mask & 0x30303030303ULL);
        v_out[2] = r2 ^ (e_mask & 0x30303030303ULL);
        v_out[3] = r3 ^ (e_mask & 0x30303030303ULL);
        // 134 lines of code, 282 operations
        v_in++;
        v_out++;
        v_in += 3;
        v_out += 3;
        }
        // End of automatically generated matrix operation.
 
    }
}  


// %%EXPORT p
void mm_op3_xi(uint_mmv_t *v_in,  uint32_t exp, uint_mmv_t *v_out)
{
    uint_mmv_t i, exp1;
 
    exp %= 3;
    if (exp == 0) {
        for (i = 0; i < 7734; ++i) v_out[i] = v_in[i];
        return;
    }
    exp1 =  (uint_mmv_t)exp - 0x1ULL;

    // Do monomial part, i.e. tags B, C, T, X
    mm_op3_xi_mon(v_in, exp1, v_out);

    // Do tag A
    mm_op3_xi_a(v_in, exp1, v_out); 

    // Do tags X, Y
    for (i = 0; i < 4; ++i) {
        uint_mmv_t *p_src = v_in + MM_OP3_OFS_Z + (i << HALF_YZ_SHIFT);
        mm_op3_xi_yz(p_src, exp1, v_out + TAB3_XI64_OFFSET[exp1][i]);
    }
    
}


