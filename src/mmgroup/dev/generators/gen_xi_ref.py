r"""Reference implementation for module ``gen_xi_functions.c``

The methods in class ``GenXi`` in this module are a reference
implementations of the functions in module ``gen_xi_functions.c``.
These functions deal with the operation of the generator
:math:`\xi` of the monster defined  in :cite:`Seysen20`,
Section 9, on a subgroup :math:`Q_{x0}`, and also with the monomial
part of that operation on the 196884-dimensional representation
:math:`\rho` of the monster.

Generator :math:`\xi` is an element of the subgroup
:math:`G_{x0}` of the monster. It operates on the normal subgroup
:math:`Q_{x0}` of :math:`G_{x0}` of structure :math:`2^{1+24}`
by conjugation as described in cite:`Seysen20`, Lemma 9.5.
Method ``gen_xi_op_xi`` of class ``GenXi`` implements the
operation of the generators :math:`\xi^e` on the group
:math:`Q_{x0}`. Here the elements of :math:`Q_{x0}` are given
in *Leech lattice encoding*, as described in 
*The C interface of the mmgroup project*, Section 
*Description of the mmgroup.generators extension*.

Representation :math:`\rho` has a subspace :math:`98280_x` of 
dimension 98280 on which :math:`G_{x0}` operates monomially in the 
same way as it operates on the short elements of :math:`Q_{x0}`. 
Here the basis vectors of  :math:`98280_x` (together with their
opposite vectors) can be identified with the short elements of
:math:`Q_{x0}`.

For a fast implementation of the operation of :math:`\xi^e` on
:math:`98280_x` we precompute tables as described in the
*mmgroup guide for developers*, Section
*Some mathematical aspects of the implementation*, Subsection
*Implementing generators of the Monster group*, Subsubsection
*Monomial operation of the generators \xi^e*.
We adopt the terminology from that subsection.

In that subsection the basis vectors of :math:`98280_x` are
subdivided into five boxes, so that :math:`\xi` acts as a
permutation on these boxes. Within each of these boxes the
basis vectors of :math:`98280_x` are numbered in a natural way
as described in that subsection. Operator :math:`\xi` permutes
the five boxes 1,...,5 as follows:

   *  1 -> 1,  2 -> 2,  3 -> 4 -> 5 -> 3 .
  
For each short vector in :math:`Q_{x0}` we encode the number of
the box containing the vector, the position of the vector inside
the box, and and the sign of the vector in the lowest 19 bits
of a 32-bit integer as follows:

  *  Bit 18...16:  number of the box, running from 1 to 5.
  
  *  Bit 15:       sign bit 
  
  *  Bit 14...0:   Entry inside the box

We call this encoding the *Short vector encoding* of a short vector
in  :math:`Q_{x0}`.

Methods ``gen_xi_short_to_leech`` and
``gen_xi_leech_to_short`` in class ``GenXi`` convert elements of
:math:`Q_{x0}` from *Leech lattice encoding* to
*Short vector encoding* and vice versa.
Method ``gen_xi_op_xi_short`` uses these two methods together with
method ``gen_xi_op_xi`` in order to compute the operation of
:math:`\xi^e` on a vector in :math:`98280_x` encoded in
*Short vector encoding*.

Using method ``gen_xi_op_xi_short`` we may easily compute tables for
the operation of :math:`\xi` and :math:`\xi^2` on a specific box.
Method ``gen_xi_make_table`` computes an array containing the lower
16 bits of the images of the entries of a box under the operation
:math:`\xi^e`. Bits 18...16 of these images can be reconstructed
from the permutation of the boxes given above.
"""
from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

from mmgroup.bitfunctions import bw24, lmap

import random
import numpy
from numpy import array, zeros, uint8, uint16, uint32

from mmgroup.dev.mat24.mat24_ref import Mat24


#######################################################################
# Auxiliary functions 
#######################################################################


def w2_gamma(v):
    """Implement function gamma() in [Seysen20], section 3.3.

    Given a Golay code vector v  in 'gcode' representation, see
    module mat24.py, we return a pair w2, c. Here c is equal to
    the element gamma(v) of the Golay cocode, with gamma() as in 
    [Seysen20], section 3.3. ,c is returned in 'cocode' representation, 
    see module mat24.py. We also return the bit w2 = w2(c), with w2()
    as defined in [Seysen20], section 3.3.   
    """
    v = int(v)
    x = Mat24.gcode_to_vect(v)
    x1 = sum( (x >> i) & 0x111111 for i in range(1,4) )
    x1 = (x1 >> 1) & 0x111111
    w2 = bw24(x1)
    w2 = ((w2 * (w2 - 1)) >> 1) & 1
    c = Mat24.vect_to_cocode(x1)
    return w2, c



def compress_gray(x):
    """Extract the gray part of a Colay code or cocode word

    The function returns the 'gray' part of a Golay code or cocode
    as a 6-bit number. It drops the 'colored' part of the word.

    A 12-bit Golay code or cocode word can be expressed as a sum of
    a 6-bit 'gray' and a 6-bit 'colored' word, see [Seysen20], 
    section 2.2.  

    A Golay code word x must be given in 'gcode' and a cocode word 
    must be given in 'cocode' representation. 
    """
    x = int(x)
    return (x & 0x0f) + ((x >> 6) & 0x30)

def expand_gray(x):
    """Reverse the effect of function compress_gray()

    Given a compressed 'gray' Golay code or cocode word, the function 
    returns the original 'gray' word in code or cocode representation.
    """
    x = int(x)
    return (x & 0x0f) + ((x & 0x30) << 6) 



#######################################################################
# Class GenXi
#######################################################################



class GenXi(object):
    """Operation of xi of a group G_x0 on the normal subgroup Q_x0.

    A subgroup G_x0 of shape 2**{1+24}.Co_1 of the Monster group is 
    used for the construction of the monster group M in [Con85]. G_x0 
    is generated by a subgroup N_x0 of shape 2**{1+24}.2**{11}.Mat24, 
    and an additional element xi of order 3, see [Seysen20], section 9. 
    Mat24 is the Mathieu group acting on 24 elements. Mat24 (and also, 
    to some extent, N_x0) is supported by class mat24.Mat24. 

    G_x0 has an extraspecial normal subgroup Q_x0 of shape 2**{1+24},
    see [Con85], [Seysen20]. Class Mat24Xi focusses on the operation of 
    xi and x**2 on Q_x0 by conjugation.

    This class provides a table-based implementation of the operation
    of xi and xi**2, which is based on [Seysen20], Lemma 9.5. The class 
    also provides some auxiliary functions used for that operation 
    defined in [Seysen20].
    
    The so-called short elements of Q_x0 (see [Con85], sect.7) play an 
    important role as basis vectors in the representation of the 
    monster in [Con85] and in this project. Class Mat24Xi also contains 
    functions for mapping short vectors in Q_x0 to basis vectors of 
    that representation and vice versa. The basis for our representa-
    tion of the monster group, is described in module mm_aux.c.

    A faster C implementation of this class is given in class
    mat24fast.Mat24Xi      

    References
    ----------
    see file refereces.txt
    """
    tab_g_gray = numpy.zeros(64, dtype = numpy.uint8)
    tab_g_cocode = numpy.zeros(64, dtype = numpy.uint8)
    #tab_g_col = numpy.zeros(64, dtype = numpy.uint8)

    tables = {
            "GenXi_module_doc" : globals()["__doc__"],
            "GenXi_g_gray_table" : tab_g_gray,
            "GenXi_g_cocode_table" : tab_g_cocode,
    }

    directives = {}

    for x in range(64):
        w2, c = w2_gamma(expand_gray(x))
        cx = compress_gray(c)
        w2x = w2 << 7 
        tab_g_gray[x] = w2x + cx 
        tab_g_cocode[cx] = w2x + x


    

    @classmethod
    def gen_xi_g_gray(cls, v):
        """Implement function gamma() in [Seysen20], section 3.3.

        Given a Golay code vector v  in 'gcode' representation, see
        module mat24.py, we return the element gamma(v) of the Golay 
        cocode, with gamma() as in [Seysen20], section 3.3. c is 
        returned in 'cocode' representation, see module mat24.py.
        """
        return expand_gray(int(cls.tab_g_gray[compress_gray(v)]))
  
    @classmethod
    def gen_xi_w2_gray(cls, v):
        """Implement function w2() in [Seysen20], section 3.3.

        Given a Golay code vector v  in 'gcode' representation, see
        module mat24.py, we return the bit w2 = w2(c), with w2() as
        defined in [Seysen20], section 3.3.   
        """
        return int(cls.tab_g_gray[compress_gray(v)]) >> 7

    @classmethod
    def gen_xi_g_cocode(cls, c):
        """Inverse of method xi_g_gray(v)

        Given a cocode vector c in in 'cocode' representation, see 
        module mat24.py, the function returns the unique gray Golay
        code vector v such that cls.gen_xi_w2_gray(v) is the gray 
        part of c. 
        """
        return expand_gray(int(cls.tab_g_cocode[compress_gray(c)]))
  
    @classmethod
    def gen_xi_w2_cocode(cls, c):
        """Implement function w2() in [Seysen20], section 3.3.

        Given a cocode vector c in in 'cocode' representation, see 
        module mat24.py, the function returns the bit w2 = w2(c), 
        with w2() as defined in [Seysen20], section 3.3.   
        """
        return int(cls.tab_g_cocode[compress_gray(c)]) >> 7

    @classmethod
    def gen_xi_op_xi(cls, x, exp):
        """Operation of  xi**exp  on the element x of the group Q_x0.

        The function returns the element

            xi**(-exp)  *  x  *  xi**exp

        of the group Q_x0. The element x of Q_x0 as must be given in
        Leech lattice encoding, as described in the header of this 
        module.

        The returned result is is coded in the same way. This 
        function implements the formula in [Seysen20], Lemma 9.5.
        """
        x = int(x)
        exp %= 3
        if (exp == 0): 
            return x
        scal = bw24((x >> 12) &  x & 0xc0f) & 1  
        x ^= scal << 24 # xor scalar product to sign

        tv =  int(cls.tab_g_gray[compress_gray(x >> 12)])
        w2v, gv = tv >> 7, expand_gray(tv)
        tc =  int(cls.tab_g_cocode[compress_gray(x)])
        w2c, gc = tc >> 7, expand_gray(tc)
        if (exp == 1):
            x &= ~0xc0f000  # kill gray code part
            x ^= w2c << 24  # xor w2(cocode) to sign
        else:
            x &= ~0xc0f     # kill gray cocode part
            x ^= w2v << 24  # xor w2(code) to sign
        x ^= gv         # xor g(code) to cocode 
        x ^= gc << 12   # xor g(cocode) to code
        return x




        
    @staticmethod
    def gen_xi_short_to_leech(x1):
        """Convert short vector to Leech lattice encoding.

        Both, Leech lattice and short vector encoding of a short vector 
        in Q_x are decribed in the header of this module. The function 
        returns the Leech lattice encoding of element x1 given in short
        vector encoding. 

        The function returns 0 for an illegal input x1. 
        """
        x1 = int(x1)
        box, sign, code = x1 >> 16, (x1 >> 15) & 1, x1 & 0x7fff
        octad = 0xffff
        if box == 1:
            if code < 1536:  # this is 2 * 24 * 32 
                gcode = code >= 768
                if gcode:
                    code -= 768
                gcode <<= 11
                i, j = code >> 5, code & 31
                cocode = Mat24.vect_to_cocode((1 << i) ^ (1 << j))
                if cocode == 0 or cocode & 0x800:
                    return 0
            elif code < 2496: # this is 2 * 24 * 32 + 15 * 64 
                octad = code - 1536
            else:
                return 0
        elif box == 2:
            if code >= 23040:  # this is  360 * 64     
                return 0   
            octad = code + 960  # this is 15 * 64
        elif box == 3:
            if code >= 24576:  # this is  384 * 64     
                return 0   
            octad = code + 24000  # this is (15 + 360) * 64
        elif box < 6:
            code += (box - 4) << 15
            cocode = Mat24.vect_to_cocode(1 << (x1 & 31))
            if cocode  == 0:
                return 0
            gcode = (code >> 5) & 0x7ff
            w =  Mat24.gcode_weight(gcode) ^ Mat24.scalar_prod(gcode, cocode) 
            gcode ^= (w & 1) << 11
        else:
            return 0  
        if octad < 48756:  # this is  759 * 64 
            cc = octad & 0x3f
            w = Mat24.suboctad_weight(cc)
            gcode = Mat24.octad_to_gcode(octad >> 6)
            gcodev = Mat24.gcode_to_vect(gcode)
            cocode = Mat24.suboctad_to_cocode(cc, octad >> 6)
            gcode ^=  w << 11
        cocode ^= Mat24.ploop_theta(gcode)
        return int((sign << 24) | (gcode << 12) | cocode)

    @staticmethod
    def gen_xi_leech_to_short(x1):   
        """Convert Leech lattice to short vector encoding.

        Both, Leech lattice and short vector encoding of a short vector 
        in Q_x are decribed in the header of this module. The function 
        returns the short vector encoding of element x1 given in Leech 
        lattice encoding. 

        The function returns 0 if the vector x1 is not short. 
        """
        x1 = int(x1)
        sign = (x1 >> 24) & 1
        x1 ^= Mat24.ploop_theta(x1 >> 12)
        gcodev = Mat24.gcode_to_vect(x1 >> 12) 
        cocodev = Mat24.cocode_syndrome(x1, 
            min(23, Mat24.lsbit24(gcodev))
        )
        w = Mat24.gcode_weight(x1 >> 12)
        if x1 & 0x800:
            if (Mat24.bw24(cocodev) > 1 or
                Mat24.scalar_prod(x1 >> 12, x1) !=  (w & 1)):
                return 0
            y = Mat24.lsbit24(cocodev)
            code = (x1 & 0x7ff000) >> 7 | y  
            box = 4 + (code >> 15)
            code &= 0x7fff
        else: 
            if w == 3:
                return 0
            elif w in [2,4]:
                code = Mat24.cocode_to_suboctad(x1, x1 >> 12, 0)
                if code >= 24000: # this is (15 + 360) * 64
                    code -= 24000
                    box = 3
                elif code >= 960: # this is 15 * 64
                    code -= 960
                    box = 2
                else:
                    code += 1536
                    box = 1
            else:
                y1 = Mat24.lsbit24(cocodev) 
                cocodev ^= 1 << y1  
                y2 = Mat24.lsbit24(cocodev) 
                if cocodev != (1 << y2) or y1 >= 24:
                    return 0
                code = 384 * (w & 2) + 32 * y2 + y1
                box = 1
        return  int((box << 16) | (sign << 15) | code)
  

    @classmethod
    def gen_xi_op_xi_short(cls, x, exp):
        """Operation of  xi**exp  on the element x of the group Q_x0.

        The function returns the element

            xi**(-exp)  *  x  *  xi**exp

        of the group Q_x0. The element x of Q_x0 as must be given in
        short vector encoding, as described in the header of this 
        module.

        The returned result is is coded in the same way. 
        """
        x = int(x)
        y = cls.gen_xi_short_to_leech(x)
        if (y == 0): return x
        y = cls.gen_xi_op_xi(y, exp)
        if (y == 0): return x
        y =  cls.gen_xi_leech_to_short(y)
        if (y == 0): return x
        return y

        
    @classmethod
    def make_table(cls, u_box, u_exp):
        assert 1 <= u_box <= 5 and 1 <= u_exp <= 2
        t_size = [0, 2496, 23040, 24576, 32768, 32768]
        a = numpy.zeros(32768, dtype = numpy.uint16) 
        length = t_size[u_box] 
        u_box <<= 16;
        for i  in range(length):
            a[i] = cls.gen_xi_op_xi_short(u_box + i, u_exp) & 0xffff;
        return a[:length]

    @staticmethod
    def invert_table(table, n_columns, len_result):
        assert len(table) & 31 == 0 and len_result & 31 == 0
        result = numpy.zeros(len_result, dtype = numpy.uint16) 
        for i, r in enumerate(table):
            if ((i & 31) < n_columns and (r & 0x7fff) < len_result):
                result[r & 0x7fff] = i | (r & 0x8000)
        return result


    @staticmethod
    def split_table(table, modulus):
        length = len(table)
        assert length & 31 == 0
        a = numpy.zeros(length >> 5, dtype = numpy.uint32) 
        for i in range(0, length, 32):
            sign = sum(((table[i+j] >> 15) & 1) << j for j in range(32))
            a[i >> 5] = sign
        table = (table & 0x7fff) % modulus
        return table, a






class Tables(GenXi):
    pass


