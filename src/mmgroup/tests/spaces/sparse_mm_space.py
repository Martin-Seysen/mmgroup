from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import sys
import os
import warnings
import copy
import numbers
import re
import numpy as np
from random import randint, randrange, sample
from functools import partial
from collections import defaultdict

from mmgroup.tests.groups.mgroup_n import MGroupNWord, StdMGroupN
from mmgroup.structures.abstract_rep_space import mod_rand_invertible
from mmgroup.structures.abstract_rep_space import mod_rand_unit
from mmgroup.structures.abstract_rep_space import mod_pwr2
from mmgroup.structures.abstract_mm_rep_space import AbstractMmRepVector
from mmgroup.structures.abstract_mm_rep_space import AbstractMmRepSpace
from mmgroup.structures.abstract_mm_rep_space import add_vector
from mmgroup.structures.mm_space_indices import purge_sparse_entry
from mmgroup.structures.parse_atoms import AtomDict
from mmgroup.structures.parse_atoms import eval_atom_expression  
from mmgroup.structures.parse_atoms import ihex 
from mmgroup.tests.spaces.rep_aux import pm_mat_from_function 
from mmgroup.tests.spaces.rep_aux import pm_diag_from_function 


from mmgroup import mat24 as m24
from mmgroup import generators as gen
from mmgroup.generators import gen_xi_w2_gray as w2_gray


######################################################################
# Modelling a vector
######################################################################

class SparseMmVector(AbstractMmRepVector):
    group = StdMGroupN
    __slots__ =  'p', 'data'
    ERR_P = "Illegal modulus %d for MM vector space"

    def __init__(self, p, tag = 0, i0 = None, i1 = None):
        if (p & 1) == 0 or not 2 < p < 256:
            raise ValueError(self.ERR_P % p)
        self.p = p
        self.data = defaultdict(int)
        add_vector(self, tag, i0, i1)

    def set_zero(self, p):
        if (p & 1) == 0 or not 2 < p < 256:
            raise ValueError(self.ERR_P % p)
        self.p = p
        self.data = defaultdict(int)

    def __len__(self):
        return len(self.data) 
  

######################################################################
# Monomial operation of the group
######################################################################

p, x, y, t, l = 0x10, 0x20, 0x30, 0x40, 0x50

gtag_dict = { 'p':p, 'x':x, 'y':y, 't':t, 'l':l}

A, B, C, T, X, Z, Y = 1, 2, 3, 4, 5, 6, 7

######################################################################
# Operation of pi (tag = 'p')

def mul_Xp(tag, d, i, g):
    eps, pi, rep =  g.cocode, g.perm, g.rep   
    d1 = m24.op_ploop_autpl(d, rep)
    i1 = pi[i]
    s = d1 >> 12
    if eps & 0x800:  # if eps is odd:
        s ^= m24.scalar_prod(d, m24.vect_to_cocode(1 << i))
        s ^= m24.pow_ploop(d, 2) >> 12
    return s, tag, d1 & 0x7ff, i1


def mul_Yp(tag, d, i, g):
    eps, pi, rep =  g.cocode, g.perm, g.rep   
    d1 = m24.op_ploop_autpl(d, rep)
    i1 = pi[i]
    s = (d1 >> 12) ^ (d1 >> 11)
    if eps & 0x800:  # if eps is odd:
        tag = Z
        s ^= d1 >> 11
    return s, tag, d1 & 0x7ff, i1

def mul_Zp(tag, d, i, g):
    eps, pi, rep =  g.cocode, g.perm, g.rep   
    d1 = m24.op_ploop_autpl(d, rep)
    i1 = pi[i]
    s = d1 >> 12
    if eps & 0x800:  # if eps is odd:
        tag = Y
        s ^= d1 >> 11
    return s, tag, d1 & 0x7ff, i1



def mul_Tp(tag, octad, sub, g):
    d = m24.octad_to_gcode(octad)
    c = m24.suboctad_to_cocode(sub, octad)
    eps, pi, rep =  g.cocode, g.perm, g.rep   
    d1 = m24.op_ploop_autpl(d, rep)
    c1 = m24.op_cocode_perm(c, pi)
    s = d1 >> 12
    s ^= (eps >> 11) & 1 & m24.suboctad_weight(sub) 
    octad1 = m24.gcode_to_octad(d1) 
    sub1 = m24.cocode_to_suboctad(c1, d1) & 0x3f
    return s, tag, octad1, sub1

def mul_Ap(tag, i, j, g):
    pi =  g.perm
    i1, j1 = pi[i], pi[j]
    i1, j1 = max(i1, j1), min(i1, j1)
    return  0, tag, i1, j1 

mul_Bp = mul_Ap

def mul_Cp(tag, i, j, g):
    eps, pi =  g.cocode, g.perm
    i1, j1 = pi[i], pi[j]
    i1, j1 = max(i1, j1), min(i1, j1)
    return  (eps >> 11) & 1, tag, i1, j1 



######################################################################
# Operation of x (tag = 'x')

def mul_Xx(tag, d, i, g):
    e = g.pl
    s = m24.scalar_prod(e, m24.vect_to_cocode(1 << i))
    s ^= m24.ploop_comm(d,e)
    return s, X, d & 0x7ff, i



def mul_Yx(tag, d, i, g):
    e = g.pl
    ed = m24.mul_ploop(m24.pow_ploop(e,3), d)
    s = m24.scalar_prod(e, m24.vect_to_cocode(1 << i))
    s = (ed >> 12) ^  (ed >> 11)  
    return s & 1, Y, ed & 0x7ff, i


def mul_Zx(tag, d, i, g):
    e = g.pl
    ed = m24.mul_ploop(m24.pow_ploop(e,3), d)
    s = m24.scalar_prod(e, m24.vect_to_cocode(1 << i))
    s = (ed >> 12)  
    return s & 1, Z, ed & 0x7ff, i

def mul_Tx(tag, octad, sub, g):
    d = m24.octad_to_gcode(octad)
    c = m24.suboctad_to_cocode(sub, octad)
    e = g.pl
    s = m24.ploop_comm(d, e)
    s ^= m24.scalar_prod(e, c)
    return s, tag, octad, sub

def mul_Ax(tag, i, j, g):
    return 0, tag, i, j

def mul_Bx(tag, i, j, g):
    c = m24.vect_to_cocode((1 << i) ^ (1 << j))
    e = g.pl
    s = m24.scalar_prod(e, c)
    return s, tag, i, j

mul_Cx = mul_Bx


######################################################################
# Operation of y (tag = 'y')


def mul_Xy(tag, d, i, g):
    e = g.pl
    de = m24.mul_ploop(d, e)
    s =  de >> 12
    return s, X, de & 0x7ff, i


def mul_Yy(tag, d, i, g):
    e = g.pl & 0xfff
    s = m24.scalar_prod(e, m24.vect_to_cocode(1 << i))
    s ^= m24.ploop_comm(d,e) 
    return s, Y, d & 0x7ff, i


def mul_Zy(tag, d, i, g):
    e = g.pl
    de = m24.mul_ploop(d, e)
    s = m24.scalar_prod(e, m24.vect_to_cocode(1 << i))
    s ^= de >> 12
    return s, Z, de & 0x7ff, i


def mul_Ty(tag, octad, sub, g):
    d = m24.octad_to_gcode(octad)
    c = m24.suboctad_to_cocode(sub, octad)
    e = g.pl
    c1 = m24.ploop_cap(d, e) ^ c
    sub1 = m24.cocode_to_suboctad(c1, d) & 0x3f
    s = m24.scalar_prod(e, c)
    return s, tag, octad, sub1


def mul_Ay(tag, i, j, g):
    c = m24.vect_to_cocode((1 << i) ^ (1 << j))
    e = g.pl
    s = m24.scalar_prod(e, c)
    if i == j: assert s == 0
    return s, tag, i, j 

def mul_By(tag, i, j, g):
    c = m24.vect_to_cocode((1 << i) ^ (1 << j))
    e = g.pl
    s = m24.scalar_prod(e, c)
    tag = B + s
    return s, tag, i, j

def mul_Cy(tag, i, j, g):
    c = m24.vect_to_cocode((1 << i) + (1 << j))
    e = g.pl
    s = m24.scalar_prod(e, c)
    tag = C - s
    return s, tag, i, j



######################################################################
# Operation of t (tag = 't')

def sign_t_XY(d, i):
    return m24.scalar_prod(d, m24.vect_to_cocode(1 << i))

def sign_t_YZ(d, i):
    return (m24.scalar_prod(d, m24.vect_to_cocode(1 << i))
            ^ (m24.pow_ploop(d, 2)) >> 12)

def sign_t_ZX(d, i):
    return m24.pow_ploop(d, 2) >> 12

def sign_t_0(d,i):
    return 0

XYZ_t_TAGS = [X, Y, Z, X, Y, Z]

XYZ_t_SIGNS = [
    (sign_t_0, sign_t_XY,sign_t_ZX),
    (sign_t_0, sign_t_YZ,sign_t_XY),
    (sign_t_0, sign_t_ZX,sign_t_YZ),
]

def mul_XYZt(index, tag, d, i, g):
    e = g.exp % 3
    return  XYZ_t_SIGNS[index][e](d,i), XYZ_t_TAGS[index + e], d, i

mul_Xt = partial(mul_XYZt, 0)
mul_Yt = partial(mul_XYZt, 1)
mul_Zt = partial(mul_XYZt, 2)



######################################################################
# Monomial operation of l on tags B, C and T
    
# The following dict maps keys (tag, i1, i2) to pairs
# (sign, (tag, i1, i2)), with (-1)**sign the sign of the monomial
# operation. The dictionary comprises all unit vectors for tags
# B, C and the unit vectors with tag T corresponding to even octads.
# These are the firt 375 octads with tag T 
dict_BCT = [None, {}, {}]


table_BCT = [None, (0, gen.make_table(3, 1)),
                 (1024, gen.make_table(3, 2)) ]



def entry_to_atom_box1(x):
    sign = (x >> 15) & 1
    x &= 0x7fff
    if x >= 1536:
        x1 = x - 1536
        return sign, (T, x1 >> 6, x1 & 63)
    elif x >= 768:
        return sign, (C, (x >> 5) - 24, x & 31)
    else:
        return sign, (B, (x >> 5), x & 31)

def entry_to_atom_box2(x):
    sign = (x >> 15) & 1
    return sign, (T, ((x >> 6) & 0x1ff) + 15, x & 63)


def unpack(a):
    return (a[0],) + a[1]

def make_dict_BCT():
    global dict_BCT
    #print("Computing dict...",  end = "")
    ar_T = [None, gen.make_table(1, 1), gen.make_table(1, 2)]
    for i in range(2496):
        t1 =  ar_T[1][i]
        if t1:
            key_ = entry_to_atom_box1(i)[1]
            dict_BCT[1][key_] = unpack(entry_to_atom_box1(t1))
            dict_BCT[2][key_] = unpack(entry_to_atom_box1(ar_T[2][i]))
    #print("length =", len(dict_BCT[1]), 375*64+2*24*24)
    ar_T = [None, gen.make_table(2, 1), gen.make_table(2, 2)]
    for i in range(360 * 64):
        key_ = entry_to_atom_box2(i)[1]
        dict_BCT[1][key_] = unpack(entry_to_atom_box2(ar_T[1][i]))
        dict_BCT[2][key_] = unpack(entry_to_atom_box2(ar_T[2][i]))


make_dict_BCT()


def op_l_BCT(tag, i, j, g):
    e = g.exp % 3
    if e == 0:
        return 0, tag, i, j
    try:
        return dict_BCT[e][tag, i, j]
    except KeyError:
        tag, o, s = tag, i, j
        assert tag == T, tag
        ofs, table = table_BCT[e]
        r = table[((o - 375) << 6) + s]
        sign = (r >> 15) & 1
        o = ofs + ((r >> 5) & 0x3ff)
        return sign, X, o, r & 0x1f


mul_Bl = mul_Cl = mul_Tl = op_l_BCT 


######################################################################
# Monomial operation of l on tag X


tag_ofs_shift = [None, None, None, (T,375,6), (X,0,5), (X,1024,5)]

def make_table_x_entry(exp, d):
     box = d + 4
     table = gen.make_table(box, exp)
     box = gen.gen_xi_op_xi_short(box << 16, exp) >> 16
     tag, ofs, shift  = tag_ofs_shift[box]
     mask = (1 << shift) - 1
     return tag, ofs, shift, mask, table


table_X = [None, [[],[]], [[],[]]]

def make_table_X():
    global table_X
    for exp in (1,2):
        for d in (0,1):
            table_X[exp][d] = make_table_x_entry(exp, d)

make_table_X()
    
    
def mul_Xl(tag, d, i, g):
    e = g.exp % 3
    if e == 0:
        return 0, tag, d, i
    tag, ofs, sh, mask, table = table_X[e][d >> 10]
    r = table[ ((d & 0x3ff) << 5) + i]
    sign = (r >> 15) & 1
    d = (r & 0x7fff) >> sh
    i = r & mask
    return sign, tag, d + ofs, i





######################################################################
# Dictionary of monomial operations

mult_monomial_dict = {}

for vtag0, v in enumerate("ABCTXZY"):
    vtag = vtag0 + 1
    for g in "pxytl":
        gtag = gtag_dict[g]
        try:
            mult_monomial_dict[vtag + gtag] = globals()["mul_" + v + g]
        except:
            pass


  

######################################################################
# Non-monomial operation of the group
######################################################################


def zero_vector(length, *data):
    """Return zero numpy vector of given length and type np.int32"""
    return np.zeros(length, dtype = np.int32)


######################################################################
# Non-onomial operation of t


class NonMonomialOp_t:
    """Helper for class NonMonomialOp for the group element t.

    Implements the non-momomial operation of the group element t**e
    on vectors with tag "ABC" and "T".
 
    General Interface see class MonomialOp_any. See class 
    NonMonomialOp for background.
    """
    __slots__ = "d_ABC", "d_T", "d_D", "exp", "p", "half", "eighth"
    mat3_e1 = np.array([[0, 2,  -2], [1, 1, 1], [1,  -1, -1]],
          dtype = np.int32)
    mat3_e2 = np.array([[0, 2,  2], [1, 1, -1], [-1,  1, -1]],
          dtype = np.int32)
    mat3 = [None, mat3_e1 , mat3_e2]
    diag64 = pm_diag_from_function(m24.suboctad_weight, 64)
    sym64 = pm_mat_from_function(m24.suboctad_scalar_prod, 64)
    mat64_e1 = sym64 @ diag64
    mat64_e2 = diag64 @ sym64
    mat64 = [None, mat64_e1, mat64_e2] 
    for m in mat64[1:]:
        assert (m @ m @ m == 512 * np.identity(64)).all()
    assert (mat64_e1 @ mat64_e2 == 64* np.identity(64)).all()
 

    def __init__(self, p, exp):
        self.d_ABC = defaultdict(partial(zero_vector,3))
        self.d_T = defaultdict(partial(zero_vector,64))
        self.d_D = defaultdict(int)
        self.exp = exp % 3
        self.p = p
        self.half = mod_pwr2(-1, p)
        self.eighth= mod_pwr2(-3, p)
        assert exp in (1,2)

    def load(self, tag, i, j, coeff): 
        coeff %= self.p
        if tag == T:
            self.d_T[i][j] = coeff 
        elif tag in (A, B, C):
            if i == j:
                assert tag == A
                self.d_D[i] = coeff
            else:
                assert i > j, (tag, (i,j))
                self.d_ABC[(i,j)][tag-1] = coeff
        else:
            raise ValueError("Bad tag %s" % tag)

    def op_ABC(self):
        m = self.mat3[self.exp]
        for (i,j), v in self.d_ABC.items():
            assert i > j, ("* ", v, (i,j))
            v = (v @ m) * self.half % self.p
            if v[0]: yield A, i, j, int(v[0])      
            if v[1]: yield B, i, j, int(v[1])      
            if v[2]: yield C, i, j, int(v[2])      
               
    def op_T(self):
        m = self.mat64[self.exp]
        for oct, v in self.d_T.items():
            v = (v @ m) * self.eighth % self.p
            for i in range(64):
                if v[i]: yield T, oct, i, int(v[i]) 

    def process(self):
        for i, value in self.d_D.items():
            if value:
                yield A, i, i, value
        yield from self.op_ABC()
        yield from self.op_T()


######################################################################
# Non-monomial operation of l



class NonMonomialOp_l:
    """Helper for class Space4096x for the group element l.

    Implements the non-momomial operation of the group element l**e
    on vectors with tag  "A", "Y", "Z".
 
    General Interface see class auto_group_space.UnprocessedVector. See 
    class Space4096x and its base class mvector.AutoGroup_RepSpace
    for background.
    """
    __slots__ = "d_YZ", "d_A", "exp", "p", "quater", "eighth"
    # Prepare 16 x 16 matrix for tag A
    # Put L4A = 2 * xi_24a with xi_24a as in [Seys19], section 9
    L4A = pm_mat_from_function(lambda i, j: i == j, 4)
    # Put L4B = xi_24b with xi_24b as in [Seys19], section 9
    L4B = pm_diag_from_function(lambda i: i > 0, 4)
    L4_POWERS = [L4A @ L4B, L4B @ L4A]
    mat16 = [None] + [np.kron(m4, m4) for m4 in L4_POWERS]

    # Prepare 64 x 64 matrix for tags Y and Z
    MSYM16 = pm_mat_from_function(lambda i, j: w2_gray(i ^ j), 16)
    MDIAG16 = -pm_diag_from_function(lambda i: w2_gray(i), 16)
    L16_POWERS = [MDIAG16 @ MSYM16, MSYM16 @ MDIAG16]
    powers_mat64 = zip(L16_POWERS, L4_POWERS)
    mat64 = [None] + [np.kron(i,j) for i,j in powers_mat64]
    cycles3 = [None, [2, 0, 1, 3], [1, 2, 0, 3]] 

    def __init__(self, p, exp):
        self.d_YZ = defaultdict(partial(zero_vector,64))
        self.d_A = defaultdict(partial(zero_vector,16))
        self.exp = exp % 3
        self.p = p
        self.quater = mod_pwr2(-2, p)
        self.eighth= mod_pwr2(-3, p)
        assert exp in (1,2)
        

    def load(self, tag, i, j, coeff): 
        coeff %= self.p
        if tag in (Y, Z):
            i += (tag == Y) << 11
            ih, il = i & 0xff0, i & 0xf
            jh, jl = j & 0x1c, j & 3
            self.d_YZ[(ih, jh)][4 * il + jl] = coeff 
            return
        if tag == A:
            ih, il = i & 0x1c, i & 3
            jh, jl = j & 0x1c, j & 3
            self.d_A[(ih,jh)][4 * il + jl] = coeff
            if ih == jh:
                self.d_A[(ih,jh)][4 * jl + il] = coeff
            return
        raise ValueError("Bad tag %s for operation of xi" % tag)

    def op_YZ(self):
        m = self.mat64[self.exp]
        cyc3 = self.cycles3[self.exp]
        for (ih, jh), v in self.d_YZ.items():
            v = (v @ m) * self.eighth % self.p
            k_hi = cyc3[ih >> 10]
            ih = ((k_hi) << 10) + (ih & 0x3ff)
            tag = Z + (ih >> 11)
            ih &= 0x7ff
            for il in range(16):
                for jl in range(4):
                    x = int(v[4 * il + jl]) 
                    if x: yield  tag, ih + il, jh + jl, x 

    def op_A(self):
        m = self.mat16[self.exp]
        for (ih, jh), v in self.d_A.items():
            v = (v @ m) * self.quater % self.p
            assert ih >= jh
            if ih > jh:
               for il in range(4):
                   for jl in range(4):
                       x = int(v[4 * il + jl]) 
                       if x: yield A, ih + il, jh + jl, x 
            else:
               for il in range(4):
                   for jl in range(il+1):
                       x = int(v[4 * il + jl]) 
                       if x: yield A, ih + il, jh + jl, x

    def process(self):
        yield from self.op_YZ()
        yield from self.op_A()

######################################################################
# Dictionary of non-monomial operations


non_monomial_dict = {
   "t": NonMonomialOp_t,
   "l": NonMonomialOp_l,
}


######################################################################
# The default instance of the monster group
######################################################################

default_monster_group = StdMGroupN

try:
    from mmgroup.mm_space import standard_mm_group
except ImportError:
    default_monster_group.target_group = None
    err = "C version of monster group support not implemented"
    warnings.warn(err, UserWarning)
else:
    default_monster_group.target_group = standard_mm_group 
    #default_monster_group.set_preimage(standard_mm_group, tuple)

######################################################################
# class SparseMmSpace
######################################################################




class SparseMmSpace(AbstractMmRepSpace):
    """Models the sparse representation 198884x of the monster group. 

    YET TO BE DOCUMENTED !!!

    """
    group = StdMGroupN
    vector_type = SparseMmVector
    space_name = "MVSp"

    def __init__(self):
        """Create a 196884-dimensional representation of the monster

        All calculations are done modulo the odd number p
        """
        pass

    #######################################################################
    # Creating vectors 
    #######################################################################


    def zero(self, p):
        return self.vector_type(p, 0)
        

    def copy_vector(self, v1):
        """Return deep copy of group element v1"""
        v2 = self.vector_type(v1.p, 0)
        v2.data.update(v1.data)
        return v2


    def __call__(self, p, tag = 0, i0 = None, i1 = None):
        return self.vector_type(p, tag, i0, i1)

    #######################################################################
    # getitem and setitem
    #######################################################################



    def getitems_sparse(self, v1, a_sparse):
        """Get items from vector v1

        Here we assert that v1 is a vector of this vector space and
        that 'a_sparse' is a one-dimensional numpy array of type
        np.uint32, containing the coordinates to be read from v1.
 
        The function must add the corresponding coordinate to each
        entry of the array 'a_sparse' to v1. All coordinates must be
        nonnegative and < 256. 

        A zero entry in the array 'a_sparse' should be ignored.
        """
        assert v1.space == self
        for i, sp in enumerate(a_sparse):
            index = purge_sparse_entry(sp) 
            if index & 0xffffff00:
                a_sparse[i] = sp  ^ (v1.data[index] & 0xff)
        return  a_sparse 

    def setitems_sparse(self, v1, a_sparse):
        """Set items from vector v1

        Here we assert that v1 is a vector of this vector space and
        that 'a_sparse' is a one-dimensional numpy array of type
        np.uint32, containing the coordinates and values to be 
        written to v1.
 
        The function must set the corresponding coordinates in v1. 
        All coordinates must be nonnegative and < 256. 

        A zero entry in the array 'a_sparse' should be ignored.
        """
        assert v1.space == self
        for sp in a_sparse:
            sp = purge_sparse_entry(sp) 
            index = int(sp & 0xffffff00)
            if index:                 
                v1.data[index] = int(sp & 0xff)
        self.reduce(v1)
        return v1

    def additems_sparse(self, v1, a_sparse):
        """Add a vector in sparse representation to vector v.

        This method takes a numpy array 'a_indices' of integers of dtype 
        numpy.uint32 containing the description of a vector v2 in sparse 
        representation. It computes 

             v  =  v + v2 .

        Here vector v is a standard vector in this space.
        """
        assert v1.space == self
        for sp in a_sparse:
            sp = purge_sparse_entry(sp) 
            index = sp & 0xffffff00
            if index:                 
                v1.data[index] += int(sp & 0xff)
        self.reduce(v1)
        return v1



    #######################################################################
    #  vector addition and scalar multiplication 
    #######################################################################
 
    def reduce(self, v1):
        to_delete = []
        data = v1.data
        for sp, value in data.items():
            value %= v1.p
            data[sp] = value
            if value == 0:
                to_delete.append(sp)
        for sp in to_delete:
            del data[sp]
        return v1

    def iadd(self, v1, v2):
        for sp, value in v2.data.items():
            v1.data[sp] += value
        return self.reduce(v1)
 
    def imul_scalar(self, v1, a):
        a = int(a % v1.p)
        for sp, value in v1.data.items():
            v1.data[sp] = value * a
        return self.reduce(v1)

    #######################################################################
    # Obtaining a vector in sparse format 
    #######################################################################
               

    def iter_as_sparse(self, v1):
        for key_, value in v1.data.items():
            value %= v1.p
            if value:
                yield key_ + value

    def as_sparse(self, v1):
        return np.fromiter(self.iter_as_sparse(v1), dtype = np.uint32)
           
    #######################################################################
    # Multiplying with a group element
    #######################################################################

    def imul_group_monomial(self, v, g_list, w):
        p = v.p
        g_tag_list = [(g, gtag_dict[g.tag]) for g in g_list]
        wd = w.data
        wd.clear()
        for vu, v_coord in v.data.items():
            sign = 0
            vtag, i, j = vu >> 25, (vu >> 14) & 0x7ff, (vu >> 8) & 0x3f
            for g, g_tag in g_tag_list:
                f = mult_monomial_dict[vtag + g_tag] 
                s, vtag, i, j = f(vtag, i, j, g)
                sign ^=  s
            wu = int((vtag << 25) + (i << 14) + (j << 8))
            wd[wu] = p - v_coord if (sign & 1) else v_coord
        
            

    def imul_group_atom(self, v, g, w):
        p = v.p
        w.data.clear()
        gtag = gtag_dict[g.tag]
        if g.tag in "tl":
            e = int(g.exp % 3)
            if e == 0:
                w.data.update(v.data)
                return
            non_monomial = non_monomial_dict[g.tag](p, e)
        for vu, v_coord in v.data.items():
            vu, v_coord  = int(vu), int(v_coord)
            vtag, i, j = vu >> 25, (vu >> 14) & 0x7ff, (vu >> 8) & 0x3f
            try:
                f = mult_monomial_dict[vtag + gtag]
            except KeyError:
                non_monomial.load(vtag, i, j, v_coord)
            else:
                sign, wtag, wi, wj = f(vtag, i, j, g)
                sign, wi, wj =  int(sign), int(wi), int(wj)  
                wu = int((wtag << 25) + (wi << 14) + (wj << 8))
                #assert 0 <= v_coord <= p                
                w.data[wu] = p - v_coord if (sign & 1) else v_coord
        if g.tag in "tl":
            for wtag, wi, wj, w_coord in non_monomial.process():
                wu = int((wtag << 25) + (wi << 14) + (wj << 8))
                w.data[wu] = w_coord
 
       
    def imul_group_word(self, v1, g):
        """Return product v1*g of vector v1 and group word g.

        v1 may be destroyed.

        This method is called for elements v1 of the space
        'self' and for elements g of the group 'self.group' only.
        """
        g_mon = []
        v, w = v1, self.zero(v1.p)
        p = v1.p
        data = v.data
        for i, coord in data.items():
            data[i] = coord % p
        for atom in v1.group(g).iter_generators():
            if not atom.tag in "tl":
                g_mon.append(atom)
            else:
                if len(g_mon):
                    self.imul_group_monomial(v, g_mon, w)
                else:
                    v, w = w, v
                g_mon = []
                self.imul_group_atom(w, atom, v)
        if len(g_mon):
            self.imul_group_monomial(v, g_mon, w)
            v, w = w, v
        return v

    #######################################################################
    # Testing equality of vectors
    #######################################################################


    def equal_vectors(self, v1, v2):
        """Return True iff vectors v1 and v2 are equal 

        This method is called for elements v1 and v2 of the space
        'self' only.
        """
        v1 = v1.reduce()
        v2 = v2.reduce()
        return v1.p == v2.p and v1.data == v2.data
        
        
StdSparseMmSpace = SparseMmSpace()
SparseMmVector.space = StdSparseMmSpace

def SparseMmV(p):
    return partial(SparseMmVector, p)

