r"""

The functions in module ``gen_xi_functions.c`` deal with the 
operation of the generator :math:`\xi` of the monster defined 
in :cite:`Seysen20`, section 9, on a subgroup :math:`Q_{x0}` and
also with the monomial part of that operation on the  
196884-dimensional representation :math:`\rho` of the monster.

Generator :math:`\xi` is an element of the subgroup 
:math:`G_{x0}` of the monster. :math:`G_{x0}` operates on its 
normal subgroup :math:`Q_{x0}` of  of structure :math:`2^{1+24}`
by conjugation as described in the **Guide**, section 
:ref:`computation-leech2`. :math:`Q_{x0}` is extraspecial 
and has center :math:`\{\pm1\}`. We have 
:math:`Q_{x0} / \{\pm 1\} \cong \Lambda / 2 \Lambda`. An element 
:math:`x` of :math:`Q_{x0}` is short if it corresponds to
a vector in the Leech lattice :math:`\Lambda` of norm 4. 
:math:`G_{x0}` also operates on the :math:`2 \cdot 98280` short 
elements of :math:`Q_{x0}`.

Representation :math:`\rho` has a subspace :math:`98280_x` of 
dimension 98280 on which :math:`G_{x0}` operates monomially in the 
same way as it operates on the short elements of :math:`Q_{x0}`. 
So the basis vectors of  :math:`98280_x` (together with their opposite 
vectors) can be identified with the short elements of :math:`Q_{x0}`.

The functions in module ``gen_xi_functions.c`` are also used for 
computing tables describing the monomial operation of 
:math:`\xi^e, e=1,2` on :math:`98280_x`.

Module ``mmgroup.dev.generators.gen_xi_ref`` is a pure python 
substitute for this set of C functions; but calculating the tables
with python takes a rather long time.

In  section :ref:`mmrep-label` we use the following names and tags 
for the basis vectors of :math:`98280_x`.


    ================ === =================================== ========
    Name             Tag Entries                             Remarks
    ================ === =================================== ========
    :math:`X^+_{ij}`  B  i, j;  0 <= j < i < 24              (1)
    :math:`X_{ij}`    C  i, j;  0 <= j < i < 24              (1)
    :math:`X_{o,s}`   T  o, s;  0 <= o < 759, 0 <= s < 64    (2)
    :math:`X_{d,j}`   X  d, j;  0 <= d < 2**11, 0 <= j < 24  (1),(3)
    ================ === =================================== ========

Remarks

(1)  i and j, 0 <= i,j < 24  refer to basis vectors of the
     Boolean vector space in which the Golay code is defined.
(2)  o is one of 759 octads, s is one of 64 even subsets of
     octad d (modulo its complement in d), as described in
     section :ref:`octads_label`.
(3)  d with 0 <= d < 2048 refers to the Parker loop element d. 
     For larger values d we have:
     
     .. math::
     
        X_{2048+k,j} = X_{k,j}, \; Y_{2048+k,j} = -Y_{k,j},  \;
        Z_{2048+k,j} = Z_{k,j},
        
        U_{4096+k,j} = -U_{k,j},  \quad  \mbox{for} \;  U = X, Y, Z;
        \; 0 <= k < 2048 \, .

We group these basis vectrors into 5 boxes (labelled 1,...,5) with 
each box containing at most :math:`3 \cdot 2^{13}` entries. Element 
:math:`\xi` permutes these boxes as follows:

   Box1 -> Box1,  Box2 -> Box2,  Box3 -> Box4 -> Box5 -> Box3 .  (1)

The mapping from the basis vectors to entries in boxes is:

   ===================== ===     =======================
   Basis vector          Box     Entry
   ===================== ===     =======================
   B[i,j]                 1         0 +  32 * i + j
   C[i,j]                 1       768 +  32 * i + j
   T[o,s], o < 15         1      1536 +  64 * o + j             
   T[o,s], 15 <= o < 375  2        64 * (o - 15) + j           
   T[o,s], o >=  375      3        64 * (o - 375) + j
   X[d,j], d <  1024      4        32 *  d + j
   X[d,j], d >= 1024      5        32 * (d - 1024) + j
   ===================== ===     =======================

This subdivision looks weird, but is has quite a few advantages:

  * The lower index (j or s) has stride 1, and the the stride of the
    higher index (i, o or d) is the lowest possible power of two. So 
    accessing an entry is easy. In the C code for a representation
    derived from :math:`\rho` the entries of a vector of such a 
    representation are strided in the same way. So the tables 
    computed in this module can be used in the C code for operator
    :math:`\xi`.

  * An entry in a box is always less than :math:`2^{15}`, so any 
    entry can be stored in a 16-bit integer, together with a sign bit.  

  * Boxes are permuted as above, so 4 tables of 16-bit integers
    and size <=  :math:`3 \cdot 2^{15}` are sufficent to encode 
    the operation of :math:`\xi`.

We remark that boxes 3 and 5 contain the short vectors with odd
Golay code words; the other boxes contain those with even
code words. Here the parity (even/odd) of a Golay code word means 
the scalar product with the cocode word :math:`\omega` which 
has six 1-bits  in column 0 of the MOG, see :cite:`Seysen20`, 
section 2.2. Short vectors with an odd cocode element occur 
precisely in boxes 4 and 5.

A similar, but finer  subdivision of the whole space :math:`\rho`
(and not only of the subspace :math:`98280_x`) is given in 
:cite:`Iva09`, section 3.4.


The operation of the generator :math:`\xi` on the group :math:`Q_{x0}`
......................................................................

Using :cite:`Seysen20`, Lemma 9.5, we can easily compute the
values :math:`\xi^{-e} x \xi^e` for all short :math:`x \in Q_{x0}`
and :math:`e = 1,2`. Here we represent an  :math:`x \in Q_{x0}`
as a 25-bit integer in Leech lattice encoding as described above.

Let :math:`\mathcal{C}` and :math:`\mathcal{C}^*` be the Golay 
code and its cocode. In :cite:`Seysen20`, section 2.2 we give 
decompositions
:math:`\mathcal{C} = \mathcal{G} \oplus \mathcal{H}` and
:math:`\mathcal{C}^* = \mathcal{G}^* \oplus \mathcal{H}^*`.

The formulas in ibid., Lemma 9.5 use auxiliary functions 
:math:`\gamma : \mathcal{G} \rightarrow \mathcal{G^*}` and
:math:`w_2 : \mathcal{G} \cup \mathcal{G^*} \rightarrow
\mathbb{F}_2` defined in ibid., section 3.3. Using the
decompositions above, we may extend :math:`\gamma` to
:math:`\mathcal{C}` and :math:`w_2` to
:math:`\mathcal{C} \cup \mathcal{C^*}`.



Short vector encoding of the short vectors in :math:`Q_{x0}`
.............................................................
  
For each short vector in :math:`Q_{x0}` we compute the number of 
the box and the entry in the box using the information given in (2). 
Then we store the vector in the lowest 19 bits of a 32-bit integer 
as follows:

  *  Bit 18...16:  number of the box
  
  *  Bit 15:       sign bit 
  
  *  Bit 14...0:   Entry in the box

Tables for the operation of :math:`\xi` and :math:`\xi^2` contain the 
lower 16 bits of that encoding. Bits 18...16 can be reconstructed 
from the permutation of the boxes given by (1).

We precomute a set of tables containing  :math:`2 \cdot 98260`
values :math:`\xi^{-e} x \xi^e`, where the short element 
:math:`x` of :math:`Q_{x0}`  is given in 
**Short vector encoding**.


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

