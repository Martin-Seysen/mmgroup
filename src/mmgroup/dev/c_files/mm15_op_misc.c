/////////////////////////////////////////////////////////////////////////////
// This C file has been created automatically. Do not edit!!!
/////////////////////////////////////////////////////////////////////////////


#include "mm_op15.h"


// %%EXPORT p
uint32_t mm_op15_copy(uint_mmv_t *mv1, uint_mmv_t *mv2)
// Copy mv1 to mv2. Here mv1 and mv2 are vectors of the
// monster group representation modulo 15.
{
    uint_fast32_t len = 15468; 
    do {
       *mv2++ = *mv1++;
    } while(--len);
    return 0; 
}



// %%EXPORT p
uint32_t mm_op15_compare(uint_mmv_t *mv1, uint_mmv_t *mv2)
//  Compare two vectors of the monster group representation.
//  Comparison is done modulo 15.
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
        tr = (((t) << 1) & 0xeeeeeeeeeeeeeeeeULL)
               | (((t) >> 3) & 0x1111111111111111ULL);
        c = t ^ tr;   // so c = 0 iff a = +- b  (mod p)
        // In case c != 0 we already know that a != b holds.
        // So assume c == 0 and hence a = +-b, i.e.  t in [0, p].
        // Then a == b (mod p) iff t == 0 or t & a in [0, p].
        // Thus is suffices to check t & a in [0, p]. 
        t &= a;
        // %%MMV_ROTL t, 1, tr
        // Put tr = t rotated left by 1 for all fields
        tr = (((t) << 1) & 0xeeeeeeeeeeeeeeeeULL)
               | (((t) >> 3) & 0x1111111111111111ULL);
        if (c | (t ^ tr)) return 1;
    } while (--len);
    return 0; 
}
   
    

// %%EXPORT p
void mm_op15_vector_add(uint_mmv_t *mv1, uint_mmv_t *mv2)
//  Vector addition: put mv1 = mv1 + mv2.
{
    uint_fast32_t len = 15468;
    uint_mmv_t a1, b1;
    uint_mmv_t a2;
    do {
        a1 = *mv1;
        b1 = *mv2++;
        a2 = ((a1 >> 4) & 0xf0f0f0f0f0f0f0fULL)
           + ((b1 >> 4) & 0xf0f0f0f0f0f0f0fULL);
        a1 = (a1 & 0xf0f0f0f0f0f0f0fULL)
           + (b1 & 0xf0f0f0f0f0f0f0fULL);
        a1 = (a1 & 0xf0f0f0f0f0f0f0fULL) 
              + ((a1 >> 4) & 0x101010101010101ULL);
        a2 = (a2 & 0xf0f0f0f0f0f0f0fULL) 
              + ((a2 >> 4) & 0x101010101010101ULL);
        a1 = a1 + (a2 << 4);
        *mv1++ = a1;
    } while (--len);
}



// %%EXPORT p
void mm_op15_scalar_mul(int32_t factor, uint_mmv_t *mv1)
//  Vector addition: put mv1 = mv1 + mv2.
{
    uint_fast32_t len = 15468;
    uint_mmv_t a1, a2;
    factor %= 15;
    if (factor < 0) factor += 15;
    do {
        a1 = *mv1;
        a2 = ((a1 >> 4) & 0xf0f0f0f0f0f0f0fULL);
        a1 = (a1 & 0xf0f0f0f0f0f0f0fULL);
        a1 *= factor;
        a1 = (a1 & 0xf0f0f0f0f0f0f0fULL) 
              + ((a1 >> 4) & 0xf0f0f0f0f0f0f0fULL);
        a1 = (a1 & 0xf0f0f0f0f0f0f0fULL) 
              + ((a1 >> 4) & 0x101010101010101ULL);
        a2 *= factor;
        a2 = (a2 & 0xf0f0f0f0f0f0f0fULL) 
              + ((a2 >> 4) & 0xf0f0f0f0f0f0f0fULL);
        a2 = (a2 & 0xf0f0f0f0f0f0f0fULL) 
              + ((a2 >> 4) & 0x101010101010101ULL);
        a1 = a1 + (a2 << 4);
        *mv1++ = a1;
    } while (--len);
}


