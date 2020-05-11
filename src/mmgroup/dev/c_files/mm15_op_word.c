/////////////////////////////////////////////////////////////////////////////
// This C file has been created automatically. Do not edit!!!
/////////////////////////////////////////////////////////////////////////////

// %%COMMENT
// TODO: comment this!!!!

#include <string.h>
#include "mm_op15.h"   
   

static uint32_t mm_op15_group_n_swap(uint_mmv_t *v, uint32_t *g, uint_mmv_t *work)
// Multiplies the vector v of the rep 106884x with the element g of 
// the subgroup N_0 of the monster group. Here g is represented as 
// an array of five 32-bit integers as in module mm_group_n.c.
// The function requires a work buffer 'work' of the same type and 
// size as the vector v. The result g*v of the group operation is 
// stored either in v or in 'work'. The fucntion returns 0 if the
// result is stored in v and 1 if it is stored in 'work'.
{
    uint_mmv_t *p0 = v, *p1 = work, *pt;
    uint_fast32_t d = g[3];
    if (g[0]) {
        mm_op15_t(p0, g[0], p1);
        pt = p0; p0 = p1; p1 = pt;
    }
    if (g[1] | g[2]) {
        mm_op15_xy(p0, g[1], g[2], g[3], p1);
        d = 0;
        pt = p0; p0 = p1; p1 = pt;
    }
    if (g[4]) {
        mm_op15_pi(p0, d, g[4], p1);
        pt = p0; p0 = p1; p1 = pt;
    } else if (d) {
        mm_op15_delta(p0, d, p1);
        pt = p0; p0 = p1; p1 = pt;
    }
    return (p0 != v);
}


// %%EXPORT p
void mm_op15_group_n(uint_mmv_t *v, uint32_t *g, uint_mmv_t *work)
{
    if (mm_op15_group_n_swap(v, g, work)) mm_op15_copy(work, v);
}


// %%EXPORT p
void mm_op15_word(uint_mmv_t *v, uint32_t *g, int32_t len_g, int32_t e, uint_mmv_t *work)
{
    uint_mmv_t *p0 = v, *p1 = work, *pt;
    uint_fast32_t round, i;
    int32_t i_start = 0, i_stop = len_g, i_step = 1;
    uint32_t sign = 0;
    uint32_t gn[5], pending;
    mm_group_n_clear(gn);
    if (e < 0) {
        i_start = len_g - 1; i_stop = i_step = -1;
        sign = 0x80000000; e = -e; 
    }
    for (round = 0; round < e; ++round) {
        for (i = i_start; i != i_stop; i += i_step) {
            pending = mm_group_n_mul_atom(gn, g[i] ^ sign);
            if (pending) {
                if (mm_op15_group_n_swap(p0, gn, p1)) {
                    pt = p0; p0 = p1; p1 = pt;
                }
                mm_group_n_clear(gn);
                if ((pending & 0x70000000) == 0x60000000) {
                    pending &= 3;
                    if (pending) {
                        mm_op15_xi(p0, pending, p1); 
                        pt = p0; p0 = p1; p1 = pt;
                    }
                } 
            }
        }
    }
    if (mm_op15_group_n_swap(p0, gn, p1)) {
       pt = p0; p0 = p1; p1 = pt;
    }
    if (p0 != v) mm_op15_copy(work, v);
}

