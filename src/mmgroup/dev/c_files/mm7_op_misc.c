/////////////////////////////////////////////////////////////////////////////
// This C file has been created automatically. Do not edit!!!
/////////////////////////////////////////////////////////////////////////////


#include "mm_op7.h"


// %%EXPORT p
uint32_t mm_op7_copy(uint_mmv_t *mv1, uint_mmv_t *mv2)
// Copy mv1 to mv2. Here mv1 and mv2 are vectors of the
// monster group representation modulo 7.
{
    uint_fast32_t len = 15468; 
    do {
       *mv2++ = *mv1++;
    } while(--len);
    return 0; 
}



// %%EXPORT p
uint32_t mm_op7_compare(uint_mmv_t *mv1, uint_mmv_t *mv2)
//  Compare two vectors of the monster group representation.
//  Comparison is done modulo 7.
//  The function returns 0 in case of equality and 1 otherwise.
{
    uint_fast32_t len = 15468;
    uint_mmv_t a, b, t, tr, c;
    do {
        a = *mv1++;
        b = *mv2++;
        // Next we compare a and b. 
        // Idea for p = 2**k-1 and unsigned k-bit integers a, b:
        // t is in [0, p] iff t == (t right rotated by one). 
        // We have a = +- b (mod p)  iff  a ^ b in [0, p].
        t = a ^ b;
        // %%MMV_ROTL t, 1, tr
        // Put tr = t rotated left by 1 for all fields
        tr = (((t) << 1) & 0x6666666666666666ULL)
               | (((t) >> 2) & 0x1111111111111111ULL);
        c = t ^ tr;   // so c = 0 iff a = +- b  (mod p)
        // In case c != 0 we already know that a != b holds.
        // So assume c == 0 and hence a = +-b, i.e.  t in [0, p].
        // Then a == b (mod p) iff t == 0 or t & a in [0, p].
        // Thus is suffices to check t & a in [0, p]. 
        t &= a;
        // %%MMV_ROTL t, 1, tr
        // Put tr = t rotated left by 1 for all fields
        tr = (((t) << 1) & 0x6666666666666666ULL)
               | (((t) >> 2) & 0x1111111111111111ULL);
        if (c | (t ^ tr)) return 1;
    } while (--len);
    return 0; 
}
   
    

// %%EXPORT p
void mm_op7_vector_add(uint_mmv_t *mv1, uint_mmv_t *mv2)
//  Vector addition: put mv1 = mv1 + mv2.
{
    uint_fast32_t len = 15468;
    uint_mmv_t a1, b1;
    do {
        a1 = *mv1;
        b1 = *mv2++;
        a1 = (a1 & 0x7777777777777777ULL) 
              + (b1 & 0x7777777777777777ULL);                     
        a1 = (a1 & 0x7777777777777777ULL) 
              + ((a1 >> 3) & 0x1111111111111111ULL);
        *mv1++ = a1;
    } while (--len);
}



// %%EXPORT p
void mm_op7_scalar_mul(int32_t factor, uint_mmv_t *mv1)
//  Vector addition: put mv1 = mv1 + mv2.
{
    uint_fast32_t len = 15468;
    uint_mmv_t a1, a2;
    factor %= 7;
    if (factor < 0) factor += 7;
    do {
        a1 = *mv1;
        a2 = ((a1 >> 4) & 0x707070707070707ULL);
        a1 = (a1 & 0x707070707070707ULL);
        a1 *= factor;
        a1 = (a1 & 0x707070707070707ULL) 
              + ((a1 >> 3) & 0x707070707070707ULL);
        a1 = (a1 & 0x707070707070707ULL) 
              + ((a1 >> 3) & 0x101010101010101ULL);
        a2 *= factor;
        a2 = (a2 & 0x707070707070707ULL) 
              + ((a2 >> 3) & 0x707070707070707ULL);
        a2 = (a2 & 0x707070707070707ULL) 
              + ((a2 >> 3) & 0x101010101010101ULL);
        a1 = a1 + (a2 << 4);
        *mv1++ = a1;
    } while (--len);
}


