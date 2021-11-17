r"""

Module ``generators`` contains the definition of the generators of the 
monster, so that they may be used in C files. It also contains support 
for the subgroups :math:`N_{0}` of structure 
:math:`2^{2+11+2\cdot11}.(\mbox{Sym}_3 \times M_{24})` and 
:math:`G_{x0}` of structure :math:`2^{1+24}.\mbox{Co}_1` of the 
monster, as described in :cite:`Con85` and :cite:`Seysen20`. Here
:math:`M{24}` is the Mathieu group acting on 24 elements, and 
:math:`\mbox{Co}_1` is the automorphism group of the 24-dimensional
Leech lattice modulo 2.

Here we fully support the computation in the subgroup :math:`N_{0}` 
based on the generators defined in this module, so that a word in the 
generators of :math:`N_{0}` can easily be reduced to a standard form. 

We also support the operation of the group :math:`G_{x0}` on the 
Leech lattice mod 2 and mod 3 (in some cases up to sign only). For a 
full support of the subgroup :math:`G_{x0}` we also have to compute 
in a Clifford group, which is implemented in module ``clifford12``.

Our set of generators of the monster group is defined in section
:ref:`mmgroup-label`. The C implementation of this set of generators
is defined in section *Header file mmgroup_generators.h*.


The intersection :math:`N_{0} \cap G_{x0}` is a group :math:`N_{x0}` 
of structure :math:`2^{1+24}.2^{11}.M_{24}`. The group :math:`M_{24}` 
(and, to some extent, also the group :math:`N_{x0}`) is 
supported by the C functions in file ``mat24_functions.c``. 

The Leech lattice and the extraspecial group :math:`Q_{x0}`
-----------------------------------------------------------

Let :math:`Q_{x0}` the normal subgroup of :math:`G_{x0}` of structure
:math:`2^{1+24}`. Then :math:`Q_{x0}` is an extraspecial 2 group and
also a normal subgroup of :math:`G_{x0}`. Let :math:`\Lambda` be the
Leech lattice. The quotient of :math:`Q_{x0}` by its center 
:math:`\{\pm1\}` is isomorphic to  :math:`\Lambda/2 \Lambda`, which
is the Leech lattice modulo 2.


For :math:`e_i \in Q_{x0}` let :math:`\tilde{e}_i` be the vector in 
the Leech lattice (mod 2) corresponding to :math:`\pm e_i`. Then
:math:`e_1^2 = (-1)^s` for 
:math:`s = \langle\tilde{e}_1, \tilde{e}_1 \rangle /2`,
where :math:`\langle.,.\rangle` is the scalar product in the
Leech lattice. For the commutator :math:`[e_1, e_2]` we have 
:math:`[e_1, e_2] = (-1)^t`, 
:math:`t = \langle \tilde{e}_1, \tilde{e}_2 \rangle`.


Leech lattice encoding of the elements of :math:`Q_{x0}`
--------------------------------------------------------

An element of  of :math:`Q_{x0}` can be written uniquely as 
a product :math:`x_d \cdot x_\delta`, 
:math:`d \in \mathcal{P} , \delta \in \mathcal{C}^*`, see
:cite:`Seysen20`, section 5. Here :math:`\mathcal{P}` is the 
Parker loop and :math:`\mathcal{C}^*` is the Golay cocode.
We encode the element 
:math:`x_d \cdot x_\delta` of :math:`Q_{x0}` as an integer 
:math:`x`  as follows:

.. math::

       x = 2^{12} \cdot d \oplus (\delta \oplus \theta(d)) \, .

Here elements of the Parker loop and elements of the cocode 
are encoded as integers as in section :ref:`parker-loop-label`
and :ref:`golay-label`.
:math:`\theta` is the cocycle given in section 
:ref:`basis-golay-label`, and ':math:`\oplus`' means bitwise
addition modulo 2. Note that a  Parker loop element is 13 bits 
long (with the most significant bit denoting the sign) and that 
a cocode element is 12 bits long.

From this representation of :math:`Q_{x0}` we obtain a 
representation of a vector in the Leech lattice modulo 2 by 
dropping sign bit, i.e. the most significant bit at position 24. 
A vector addition in the Leech lattice modulo 2 can be done by 
applying the XOR operator ``^`` to the integers representing 
the vectors, ignoring the sign bit.
 
Special elements of the group :math:`Q_{x0}`
--------------------------------------------

We write :math:`\Omega` for the positive element of the Parker 
loop such that :math:`\tilde{\Omega}` is the Golay code word
:math:`(1,\ldots,1)` as in :cite:`Con85` and :cite:`Seysen20`.
In this specifiction we also write :math:`\Omega` for the 
element :math:`x_{\Omega}` of :math:`Q_{x0}` and for the element 
:math:`\tilde{x}_{\Omega}` of the Leech lattice modulo 2 if the 
domain of :math:`\Omega` is clear from the context. Then
:math:`\Omega` has Leech lattice encoding ``0x800000`` in
our chosen basis of the Golay code; and the element :math:`\Omega` 
of  :math:`\Lambda/2 \Lambda` corresponds to the standard
coordinate frame of the real Leech lattice.

For fast computations in the monster group it is vital to compute
in the centralizer of a certain short element :math:`x_{\beta}`
of :math:`Q_{x0}`, where :math:`\beta` is an even coloured element 
of the Golay cocode, as described in :cite:`Seysen20`. Here we 
choose the cocode element :math:`\beta` corresponding to the
element :math:`(0,0,1,1,0,\ldots,0)` of  :math:`\mbox{GF}_2^{24}`.
Then the centralizer of :math:`x_{\beta}` is isomorphic to the
a double cover baby monster group and contains the generators 
:math:`\tau` and  :math:`\xi` of the monster. We also write  
:math:`\beta` for the element :math:`x_{\beta}` of :math:`Q_{x0}` 
and for the element :math:`\tilde{x}_{\beta}` of 
:math:`\Lambda/2 \Lambda` in the same way as for the element 
:math:`\Omega`.  Then :math:`\beta` has Leech lattice 
encoding ``0x200``  in our chosen basis.

Module ``gen_leech_reduce.c`` contains functions for rotating 
arbitrary type-4 vectors in :math:`\Lambda/2 \Lambda` to 
:math:`\Omega` and for  rotating arbitrary type-2 vectors in 
:math:`\Lambda/2 \Lambda` to :math:`\beta`. 

Computations in the Leech lattice modulo 3
------------------------------------------

For the construction of the subgroup :math:`G_{x0}` of the monster
we also require the automorphism group :math:`\mbox{Co}_0` of the
**real** Leech lattice, as decribed in :cite:`Con85`. That group has 
a faithful representation as an automophism group of 
:math:`\Lambda/3 \Lambda`, but not of :math:`\Lambda/2 \Lambda`. 
Module ``gen_leech3.c`` provides functions for computing in the 
Leech Lattice modulo 3.

Note that in :cite:`Seysen20` the operation of the generators 
:math:`x_d, y_d, x_\delta, x_\pi, \xi` of :math:`G_{x0}` is
also defined on the real Leech lattice.
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
    return (x & 0x0f) + ((x >> 6) & 0x30)

def expand_gray(x):
    """Reverse the effect of function compress_gray()

    Given a compressed 'gray' Golay code or cocode word, the function 
    returns the original 'gray' word in code or cocode representation.
    """
    return (x & 0x0f) + ((x & 0x30) << 6) 



#######################################################################
# Class Mat24Xi
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
        return expand_gray(cls.tab_g_gray[compress_gray(v)])
  
    @classmethod
    def gen_xi_w2_gray(cls, v):
        """Implement function w2() in [Seysen20], section 3.3.

        Given a Golay code vector v  in 'gcode' representation, see
        module mat24.py, we return the bit w2 = w2(c), with w2() as
        defined in [Seysen20], section 3.3.   
        """
        return cls.tab_g_gray[compress_gray(v)] >> 7

    @classmethod
    def gen_xi_g_cocode(cls, c):
        """Inverse of method xi_g_gray(v)

        Given a cocode vector c in in 'cocode' representation, see 
        module mat24.py, the function returns the unique gray Golay
        code vector v such that cls.gen_xi_w2_gray(v) is the gray 
        part of c. 
        """
        return expand_gray(cls.tab_g_cocode[compress_gray(c)])
  
    @classmethod
    def gen_xi_w2_cocode(cls, c):
        """Implement function w2() in [Seysen20], section 3.3.

        Given a cocode vector c in in 'cocode' representation, see 
        module mat24.py, the function returns the bit w2 = w2(c), 
        with w2() as defined in [Seysen20], section 3.3.   
        """
        return cls.tab_g_cocode[compress_gray(c)] >> 7

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
        exp %= 3
        if (exp == 0): 
            return x
        scal = bw24((x >> 12) &  x & 0xc0f) & 1  
        x ^= scal << 24 # xor scalar product to sign

        tv =  cls.tab_g_gray[compress_gray(x >> 12)] 
        w2v, gv = tv >> 7, expand_gray(tv)
        tc =  cls.tab_g_cocode[compress_gray(x)] 
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
            w = Mat24.bw24(cc) 
            cc = (cc << 1) ^ (w & 1)
            w += w & 1
            gcode = Mat24.octad_to_gcode(octad >> 6)
            gcodev = Mat24.gcode_to_vect(gcode)
            cocode = Mat24.vect_to_cocode(Mat24.spread_b24(cc, gcodev))
            gcode ^=  ((w >> 1) & 1) << 11
        cocode ^= Mat24.ploop_theta(gcode)
        return (sign << 24) | (gcode << 12) | cocode

    @staticmethod
    def gen_xi_leech_to_short(x1):   
        """Convert Leech lattice to short vector encoding.

        Both, Leech lattice and short vector encoding of a short vector 
        in Q_x are decribed in the header of this module. The function 
        returns the short vector encoding of element x1 given in Leech 
        lattice encoding. 

        The function returns 0 if the vector x1 is not short. 
        """   
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
                if w == 4:
                    gcodev ^= 0xffffff
                    x1 ^= 0x800000
                w_bad = (Mat24.bw24(cocodev) ^ 2 ^ w) & 3
                if w_bad or (cocodev & gcodev) != cocodev:
                    return 0
                c =  Mat24.extract_b24(cocodev, gcodev)
                if (c & 0x80):
                     c ^=  0xff 
                y1 = Mat24.gcode_to_octad(x1 >> 12)
                code = (y1 << 6) | (c >> 1)
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
        return  (box << 16) | (sign << 15) | code
  

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

