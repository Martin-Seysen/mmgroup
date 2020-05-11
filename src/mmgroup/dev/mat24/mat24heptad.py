from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import types
import sys
import re
import os
from operator import __or__, __xor__, __and__
from functools import reduce

from mmgroup.bitfunctions import v2, lmap, lrange




from mmgroup.dev.generate_c.generate_functions  import prepend_blanks
from mmgroup.dev.mat24.mat24aux import generate_golay24_decode


def lsbit24(x):
    return v2(x | 0x1000000)



class HeptadCompleter(object):
    """

    Member functions tables(), directives() and bytelen() and also
    compute() and generate() follow the conventions in 
    class lsbit24_function. 
    """
    def __init__(self, gc):
        self.gc = gc 
        self.make_tables()
        self.make_exported_tables()
        self.separator = self.separate_6_7()
        self.std_octad = gc.vect_to_octad(0xff)
        self.std_octad_table_name = "MAT24_STD_OCTAD"


    def vintern_syndrome(self, v):
        """Return Golay syndrome of vector v in internal representation        

        v must be a vector with odd parity in internal representation.
        """
        syn = self.gc.syndrome_table[(v >> 1) & 0x7ff]
        return ((1 << (syn & 31)) ^ (1 << ((syn >> 5) & 31))
                                  ^ (1 << ((syn >> 10) & 31)))


    def syndrome(self, v):
        """Return Golay syndrome of vector v in standard representation        

        v must be a vector with odd parity in standard representation.
        """
        return self.vintern_syndrome(self.gc.vect_to_vintern(int(v)))

    def syndrome_with_error(self, v):
        """Same as method syndrome, but return pair (syndrome, error)
   
        syndrome is the Golay code syndrome of the vector v given in
        standard representation.        

        error is set if v does not have odd parity.
        """
        v1 = self.gc.vect_to_vintern(int(v))
        return self.vintern_syndrome(v1), (v1 & 1) ^ 1



    def make_tables(self):
        self.start = [0,1,2,3,4,5,8]
        s = reduce(__or__, [1 << i for i in self.start])

        self.syndromes = []
        self.index_table = []
        for i in range(6):
            for j in range(i):
                self.index_table.append((i,j))
                self.syndromes.append(self.syndrome(s ^ (1 << i) ^ (1 << j)))
                
        self.find_table = [None]*24
        self.imax = 0
        for i1, s1 in enumerate(self.syndromes):
            for i2, s2 in enumerate(self.syndromes[:i1]):
                a = s1 & s2
                b = lsbit24(a)                    
                if a == 1 << b and self.find_table[b] == None:
                    self.find_table[b] = [i1, i2]
                    self.imax = max(i1+1, self.imax)
        assert not None in self.find_table[9:]
        assert self.imax <= 10  # An a posteriori finding after programming!!!
        self.index_table = self.index_table[:self.imax]


    def make_exported_tables(self):
        self.exported_index_table = [ 
           (i, j1, j2) for i, (j1,j2) in enumerate(self.index_table)
        ]
        self.exported_find_table = [ 
           (i + 9, j1, j2) for i, (j1,j2) in enumerate(self.find_table[9:])
        ]


    def separate_6_7(self):
        for i1, s1 in enumerate(self.syndromes):
            for j1 in range(6):
                for j2 in range(9,24):
                    s = (1 << j1) | (1 << j2) |  s1
                    if lsbit24(s) == 5 and  self.syndrome(s) & 0xc0 == 0x40:
                         self.imax = max(i1+1, self.imax)
                         return (j1, j2, i1)


    def compute(self, p):
        """Complete a permutation given by p to an element of  Mat24.

        p must be a list of length 24. Entries p[i], i = 0,1,2,3,4,5,8
        must make up a valid umbral heptad, i.e. a heptad not contained 
        in an octad. p[0],...,p[5] must be contained in an octad, p[8]
        must not be contained in that octad. The other entries of 
        input p are ignored.

        It can be shown that such a permutation p can be completed to 
        a unique element of Mat24.

        The function returns 0 in case of success and a nonzero value
        otherwise. In case of success, p is completed to an element of
        the Mathieu group Mat24. 
        """
        start = 1 << p[8]
        for i in range(6):
            start |= 1 << p[i]
        a = []
        for i, (j1, j2) in enumerate(self.index_table):
            x = start ^ (1 << p[j1]) ^  (1 << p[j2])
            x = self.syndrome(x)
            a.append(x)

        for i, (j1, j2) in enumerate(self.find_table[9:]):
            p[i+9] = lsbit24(int(a[j1] & a[j2]))  

        syn67, err  = self.syndrome_with_error(start)
        syn67  ^= (1 << p[8]) 
        j1, j2, i = self.separator
        s = self.syndrome((1 << p[j1]) ^ (1 << p[j2]) ^ a[i])
        p[6] = lsbit24(int(syn67 & s))
        p[7] = lsbit24(int(syn67 & ~s))

        err |= (p[0] | p[1] | p[2] | p[3] | p[4] | p[5] | p[8]) & -32
        err |= start & (0xff000000 | (1 << p[6]) | (1 << p[7]) | syn67)
        err |= (syn67 - 1) & 0xff000000
        return err
        





    def int_to_perm(self, k):
        """Return the k-th permutation in the Mathieu group Mat24   

        Any integer 0 <= k < 244823040 is evaluated in mixed-radix with bases
        759, 8, 7, 6, 5, 4, 3, 16, with valence decreasing from left to right. 
        The first digit determines the number of the octad which is the image 
        of the standard octad (0,...,7). The following six digits determine 
        an even permutation of that octad, up to a permutation of the last
        two entries. The final digit determines the image of element 8, i.e. 
        the first element not in the standard octad. 

        The images of the remaining elements 6, 7 and 9,...,23  are determined 
        by calling function self.compute()
        """
        oct, k = divmod(k, 322560)
        if oct >= 759: return None
        oct -= 759 - self.std_octad
        oct += (oct >> 12) & 759 # give number 0 to standard octad
        #print("i2p", oct)
        oct = self.gc.octad_to_vect(oct);
        #print("i2p oct", hex(oct))
        p = [None]*24
        oct, j = 8 * oct, 0x8        
        for i in range(24):
            o = oct & 8
            p[(j >> o) & 0x1f] = i
            j += 1 << o
            oct >>= 1
        p[8] = p[8 + (k & 15)]
        #print("i2pfinal", k & 15)
        k >>= 4
        k *= (1 << 28) // 2520 + 1
        for i in range(6):
            k1 = i + (k >> 28)
            #print("i2p%d" % i, k >> 28)
            p[i], p[k1] = p[k1], p[i]  
            k = (k & 0xfffffff) * (7-i)
        self.compute(p)
        return p

    def perm_to_int(self, p):
        """Convert the permutation p in the Mathieu group Mat24 to an integer.

        This reverses member function int_to_perm(). The input permutation
        is not checked.
        """
        oct = sum(1 << x for x in p[:8])
        #print("p2i oct", hex(oct))
        try:
            res = self.gc.vect_to_octad(oct) 
        except:
            return -1
        res -=  self.std_octad
        res += (res >> 12) & 759 
        #print("p2i", res)
        p1 = [24]*32
        oct, j = 8 * oct, 0x00        
        for i in range(24):
            o = oct & 8
            p1[i] = (j >> o) & 0x1f
            j += 1 << o
            oct >>= 1
        q, q_inv = [None]*8, [None]*8
        for i in range(8):
            j = p1[p[i] & 0x1f] & 7
            q[j] = i
            q_inv[i] = j
        for i in range(6):
            # exchange place i with place q_inv[i]
            j = q_inv[i]
            #q_inv[q[i]], q_inv[q[j]] = q_inv[q[j]], q_inv[q[i]]
            #q[i], q[j] = q[j], q[i]
            #assert q[:i] == q_inv[:i] == lrange(i)
            q_inv[q[i]] = q_inv[q[j]]
            q[j] = q[i]
            #print("p2i%d" % i, j-i)           
            res = res * (8 - i) + j - i
        #print("p2ifinal", p1[p[8] & 0x1f])  
        return 16 * res + p1[p[8] & 0x1f]
            
            
    def perm_from_heptads(self, h1, h2):
        """Try to find a permutation p that maps heptad h1 to h2

        h1 and h2 must be lists of length 7 defining two umbral heptads,
        i.e. heptads not contained in an octad. If a permutation p in
        the Mathieu group Mat24 that maps h1 to h2 exists, it is unique. 

        Return permutation p if such a p exists an is unique,
        and return None otherwise.
        """
        # First find the special element of v h1 not contained in the octad
        v = 0
        for i in range(7):
            v |= 1 << (h1[i] & 31)
        y = self.syndrome(v)
        v =  lsbit24(v & y)
        
        # Find position y of element v in h1
        y = 0
        for i in range(7): 
            y |= ((h1[i] != v) - 1) & i

        # Copy special element of h1 to position 8 of p1 and copy the other
        # elements of h1 to positions 0,...,6. Copy h2 similarly to p2
        p1 = h1[:7] + [None]*17
        p2 = h2[:7] + [None]*17
        p = [None] * 24
        p1[8] = p1[y]
        p1[y] = p1[6]
        p2[8] = p2[y]
        p2[y] = p2[6]

        # Complete p1 and p2 from heptad. Return error if any completion fails
        if self.compute(p1) | self.compute(p2):
            return None

        # If success, return p1**(-1) * p2
        for i in range(24):
            p[p1[i]] = p2[i]
        return p


    def tables(self):
        """Return a single table for internal use only, plus
        two tables with constants to be defined with '#define' in
        the generated code

        Entries of the table are integers of at least 8 bit length.
        """
        return {
            "Mat24_heptad_find_table": self.exported_find_table,
            "Mat24_heptad_index_table": self.exported_index_table,
            "Mat24_heptad_const": list(self.separator),
            "MAT24_STD_OCTAD" : [ self.std_octad ],
        }


    def directives(self):
        return { }

    def bytelen(self):
        return len(self.find_table)



