"""A fast random generator for the representation of the monster

A fast random generator for the representation of the monster
group is provided in C. This generator supports very fast
generation of long random vectors of integers modulo p for
p <= 256.

The seed for such a random generator is given by an array of type
uint8_t seed[MM_RNG_SIZE]. In this version we have 
MM_RNG_SIZE = 266

The function mm_rng_seed(uint8_t *seed, ...) fills an existing
array 'seed' of that type with seed data. See description of
that function for details. 

The function
mm_rng_gen_modp(uint8_t p, uint8_t *seed, uint8_t *out, uint32_t len)

writes len uniform random numbers x_i with 0 <= i < p  to the array 
referred by the pointer out. Here 1 <= p < 256 must hold.

This function is optimized for generating large random vectors
with, say, len >= 1000, as required for the representation of
the monster. 

Internal operation.

We use the RC4 random generator which generates bytes in a fast
and secure manner. The RC4 random generator also initializes a
63-bit LFSR. For generating a number moduolo p, a byte is taken
from the RC4 generator and it is mixed with 48 bits taken from
the LFSR so that we obtain an uniform distributed random
number between 0 and 1. That number is multiplied by p and
the integral part of that product is returned. The LFSR is
shifted by 32 bytes in each step.

"""

from mmgroup.dev.mm_basics.mm_basics import MM_Basics

# Export the doc string to the code generator in 

class MM_Random_Doc(MM_Basics):
    def __init__(self):
        super(MM_Random_Doc, self).__init__()
        self.make_tables( {"MM_RANDOM_DOC" : __doc__ } )
