r"""Investigate certain orbits of 2A axes of the monster group

The functions in this module are used to find the orbits of 2A axes
of the subgroup :math:`G_{x0}` (of structure 
:math:`2^{1+24}.\mbox{Co}_1`) of the monster group. According to
|Nor98| there are 12 such orbits.

"""

import sys
import os
import time
from math import floor, ceil
from random import randint, shuffle, sample
from collections import defaultdict, OrderedDict
from collections.abc import Iterable
import numpy as np
from operator import __or__
from functools import reduce

if __name__ == "__main__":
    sys.path.append("../../../")

from mmgroup import MM0, MMV, MMVector, Cocode, XLeech2, Parity, PLoop
from mmgroup.mm_crt_space import MMVectorCRT 
from mmgroup import mat24, GcVector, AutPL, Parity, GCode
from mmgroup.clifford12 import bitmatrix64_solve_equation


from mmgroup.tests.test_axes.get_sample_axes import G_CENTRAL
from mmgroup.tests.test_axes.get_sample_axes import G_AXIS
from mmgroup.tests.test_axes.get_sample_axes import V_AXIS
from mmgroup.tests.test_axes.get_sample_axes import G_AXIS_OPP
from mmgroup.tests.test_axes.get_sample_axes import V_AXIS_OPP
from mmgroup.tests.test_axes.get_sample_axes import g_central
from mmgroup.tests.test_axes.get_sample_axes import g_axis
from mmgroup.tests.test_axes.get_sample_axes import g_axis_opp

V15 = MMV(15)



v_axis = MMVectorCRT(16, V_AXIS)
v_axis_opp = MMVectorCRT(16, V_AXIS_OPP)
v_axis15 = V15(V_AXIS)
v_axis_opp15 = V15(V_AXIS_OPP)


#################################################################
# Transitive closure and identifying block of matrix A
#################################################################

def transitive_closure(A):
    """Return transitive closure of symmetric square matrix ``A``

    The function returns a 0-1 matrix ``c`` of the same size as
    matrix ``A`` if ``c[i,j] = 1`` iff there is a path of nonzero
    entries 

    ``A[k[0],k[1]], A[k[1],k[2]], ... , A[k[l-1],k[k]],``

    with ``k[0] = i, k[l] = j``. Diagoanl entrie of ``c`` are set
    to 1.
    """
    A = np.array(A != 0, dtype = np.uint32)
    assert len(A.shape) == 2 and A.shape[0] ==  A.shape[1]
    A_old = np.eye(A.shape[0], dtype = np.uint32)
    A = A | A_old
    while (A != A_old).any():
        A_old = A
        A = (A @ A) != 0
    return np.array(A, dtype = np.uint32)

def bit_count(x):
    return bin(x).count("1")

def block_sort_key(x):
    return -bit_count(x), x & -x, x 

def blocks(A):
    """Return the list of blocks of a symmetric square matrix ``A``.

    Here a block a set of indices of matrix ``A`` which are in the
    same equivalence class in the transitive closure of ``A``.

    Each block in the returned list is encoded as an integer with 
    bit ``i`` of the integer set iff index ``i`` is in that block.

    Blocks are ordered by size with larger block occuring first.
    """
    A = np.array(A != 0, dtype = np.uint32)
    A = np.array(A | A.T, dtype = np.uint32)
    A = transitive_closure(A)
    s = set()
    for row in A:
        d = int(sum(x << i for i, x in enumerate(row)))
        s.add(d)
    data = list(s)
    return sorted(data, key = block_sort_key)


def block_to_list(bl):
    """Convert block returned by function blocks() to list of integers"""
    return [x for x in range(24) if bl & (1 << x)]


def adjacent_blocks(A):
    """Return adjacent blocks of a symmetric square matrix ``A``.

    If matrix ``A`` may be split into ``n`` adjacent blocks then  
    the function returns a list ``i`` of lenggth ``n+1`` such that
    all nonzero entries of ``A`` occur in one of the submatrices

    ``A[i[0]:i[1], i[0]:i[1]], ...  A[i[n]:i[n+1], i[n]:i[n+1]]``

    Then ``i[0] = 0`` and `i[n+1]`` is equal to the size of the
    square matrix ``A``.  
    """
    bl = blocks(A)
    adjacent_blocks = []
    while len(bl):
        x = max(bl)
        adjacent_blocks.append(x.bit_length())
        xmin = x & -x
        for i, y in enumerate(bl):
            if y >= xmin:
                 del bl[i]
    adjacent_blocks.append(0)
    return list(reversed(adjacent_blocks))
    

#################################################################
# Eigenvalues of matrix A
#################################################################


def purge_diag(diag, power = 1):
    """Convert list of real numbers to a dictionary

    Given a list 'diag' of rea numbers, the function returns a
    dictionary mapping each number occuring in that list to it
    multiplicity. The function makes reasonable assumptions
    regarding equality of approximately equal entries, so that
    we obtain a dictionary from a list of eigenvalues of
    a 24 times 24 matrix acting on the Leech lattice.
    """
    EPS = 1.0e-8
    data = {}
    non_ints = []
    for x0 in diag:
        x = x0**power 
        if abs(x - round(x)) < EPS:
            i = int(round(x))
            if i in data:
                data[i] += 1
            else:
                data[i] = 1
        else: 
            done = False           
            for d in non_ints:
                if abs(d - x) < EPS:
                    data[d] += 1
                    done = True
            if not done:
                data[x] = 1
                non_ints.append(x)
    return data


def format_eigen_values(A):
   """Return description of eigenvalues of a matrix ``A``

   The description of the eigenvalues is returned as a string.
   """
   eigen =  purge_diag(np.linalg.eigvals(A))
   eigenv = sorted(eigen.items(), key = lambda x: x[0])
   data = [("%d^%d" if isinstance(k, int) else "%.2f^%d") % 
       (k,v) for k, v in eigenv ]
   return "(" + ", ".join(data) + ")"
 



def block_eigenvalues(A):
    """Return description of eigenvalues of a matrix ``A``

    The description of the eigenvalues is returned as a string.
    This is similar to the function ``format_eigen_values``. But
    if matrix ``A`` consists of blocks then the description takes
    care of that block structure.
    """
    ad = adjacent_blocks(A)
    if len(ad) > 24:
        return format_eigen_values(A)
    eigen = []
    exp = []
    for i in range(len(ad) - 1):
        m = A[ad[i]:ad[i+1], ad[i]:ad[i+1]] 
        str_eigen = format_eigen_values(m)
        if len(eigen) and eigen[-1] == str_eigen:
            exp[-1] += 1
        else:
            eigen.append(str_eigen)
            exp.append(1)
    for i, s in enumerate(eigen):
        if exp[i] > 1:
            eigen[i] += "^" + str(exp[i])
    return  ", ".join(eigen)
        

#################################################################
# Dealing with one large block (of size >= 5)
#################################################################


def beautify_single_block(A):
    """Beautify matrix ``A`` acting on the Leech lattice using diagonal

    Here ``A`` is an integer 24 times 2 times matrix acting on the 
    Leech lattice. The function tries to find a permutation ``g`` in
    the group ``M_24`` that beautifies matrix ``A``. The search for
    ``g`` is based on the diagonal elements of ``A``. Here we try
    to group equal diagonal elements together.

    Element ``g`` is returned and an instance of class ``MM``
    corresponding to an element of the monster group. The function
    returns the neutral element of the monster if beautification is
    not possible.
    """
    d = defaultdict(list)
    for i in range(24):
        d[A[i,i]].append(i)
    lists = sorted(d.values(), key = lambda x: len(x))
    l = lists[0]
    if 1 <= len(l) <= 5:
        pi = AutPL(0, zip(l, range(len(l))), 0) 
        return MM0('p', pi)
    return MM0()
      



def beautify_large_block(A):
    """Beautify block matrix ``A`` acting on the Leech lattice

    The function tries to find a permutation ``g`` in the group 
    ``M_24`` that beautifies a 24 times 24 matrix ``A`` acting on 
    the Leech lattice. The objective of that beautification is to 
    bring the lines and columns of ``A`` coressponding to the
    largest block in the matrix together.

    Here a block is a collection of sets of indices. Two indices 
    ``i, j`` belong to the same block if ``A[i,j] != 0``. The
    function processes the largest block, provided that this
    block has size at least 3.

    Element ``g`` is returned and an instance of class ``MM``
    corresponding to an element of the monster group. The function
    returns the neutral element of the monster if beautification is
    not possible.
    """
    bl = blocks(A)[0]
    v = GcVector(bl)
    if len(v) < 3:
        return MM0()
    if len(v) <= 5:
        b = v.bit_list
        pi = AutPL(0, zip(b, range(len(v))), 0)
        return MM0('p', pi)
    try:
        syn = v.syndrome()
    except:
        return MM0()
    single = (v & syn).bit_list
    v &= ~syn
    if len(v) > 8:
        return MM0()
    v_list = v.bit_list
    if len(single):
        r_list = [7-x for x in range(len(v_list))]
    else:
        r_list = list(range(len(v_list)))
    if len(single):
        v_list += single
        r_list.append(8)
    for i in range(2000):
        shuffle(v_list)
        try:
            pi = AutPL(0, zip(v_list, r_list), 0)
            return MM0('p', pi)
        except:
            pass 
    return MM0()
   







def solve_gcode_diag(l):
    """Solve cocode equation

    Here ``l`` is a list of tupeles ``(i0, i1, k)``.  For an
    unknown Golay code word ``x``, each tuple means an equation
    ``<x, Cocode([i0,i1])> = k``, where ``<.,.>`` is the scalar
    product. If a solution ``x`` exists then the function
    returns ``x`` as an instance of class |PLoop|.
    """ 
    a = np.zeros(len(l), dtype = np.uint64)
    for i, (i0, i1, k) in enumerate(l):
        a[i] = Cocode([i0,i1]).ord + ((int(k) & 1) << 12)
    v = bitmatrix64_solve_equation(a, len(l), 12)
    if v < 0:
        err = "Off-diagonal matrix equation has no solution"
        raise ValueError(err)    
    result = PLoop(v)
    for i0, i1, k in l:
        c = Cocode([i0,i1])
        check = hex(result.ord), hex(c.ord), k
        assert result & c == Parity(int(k)), check
    return result


def change_signs_A(A, plus_list = [], minus_list = [], ignore = False):
    """Adjust signs of  matrix ``A`` acting on the Leech lattice 

    Here ``plus_list`` is a pair of indices ``(i,j)`` so that we
    want to achieve ``A[i,j] >= 0``, and ``minus_list`` is a pair 
    of indices ``(i,j)`` so that we want to achieve  ``A[i,j] <= 0``. 
    Any of the entries ``i,j`` of the pair may also be a list of 
    integers.

    The function tries to compute an integer ``d`` encoding an 
    element of of the Golay code, so that the transfromation
    with that Golay code element changes the signs appropriately.

    If this fails then the functoion raises ValueError. The user
    may set parameter ``ingore`` to ``True``; then the function
    returns the integer 0 in case of failure.
    """
    eqn_list = []
    for sgn, l in zip([1,-1], [plus_list, minus_list]):
        for il, jl in l:
            if not isinstance(il, Iterable): il = [il]
            if not isinstance(jl, Iterable): jl = [jl]
            for i in il:
                for j in jl:
                    if A[i,j] and i != j:
                        eqn_list.append((i, j, sgn*A[i,j] < 0))
    try:
        return solve_gcode_diag(eqn_list)
    except ValueError:
        if (ignore): return 0
        raise
             


def try_make_blocks_positive(A, block = 0, first_row_only = False, tag = 'y'):
    """Make block of matrix ``A`` acting on the Leech lattice positive

    Here ``block`` is a block of rows of columns of the matrix 
    given as a set of integers. Block may also be a nonzero
    integer encoding a block in the bits being set in that integer.
    In case ``block = 0`` (default) we take the largest block,
    where a block is as in function ``blocks()``.

    The function tries to compute an integer ``d`` encoding an 
    element of of the Golay code, so that the transfromation
    with that Golay code element makes all entries of that 
    block nonnegative. In case of success the function returns 
    the element ``return MM0('y', d)`` of the Monster group that 
    changes the signs appropriately.

    If this fails then the functoion raises ValueError. 
    """
    if isinstance(block, Iterable):
        bl_list = list(block)
    else:
        bl = block if block else blocks(A)[0]         
        bl_list = [x for x in range(24) if (1 << x) & bl] 
    row_list = bl_list[:1] if first_row_only else bl_list
    d = change_signs_A(A, [(bl_list, row_list)])
    return MM0(tag, d)


def make_blocks_positive(A, tag = 'y'):
    """Try to make matrix ``A`` acting on the Leech lattice positive

    The function tries to make the entries of the matrix  ``A``
    acting on the Leech lattice positive. It first tries to make
    all entries nonnegative and then it tries to make the entries
    of the largest block of ``A`` nonnegative.

    The function calls function  ``try_make_blocks_positive``.
    and returns an element ``MM0('y', d)`` as returned by that
    function in case of success. It returns the neutral element
    of the monster in case of failure.
    """
    try:
        return try_make_blocks_positive(A, 0xffffff, False, tag)
    except:
        try:
            return try_make_blocks_positive(A, 0, False, tag)              
        except:
            try:
                return try_make_blocks_positive(A, 0, True, tag)              
            except:
                return MM0()





#################################################################
# Dealing with one largest block of size 3
#################################################################



def sort_single_size3_block(A):
    """Sort entries of matrix with single 3 times 3 block

    If matrix ``A`` contains a single 3 times 3 block and some
    diagonal entries the the function returns a permutation that 
    moves the 3 times 3 block to rows and columns 0, 1, and 2.
    That permutation is returned as an instance of class ``MM0``.
    
    Otherwise the function returns the neutral element in
    class ``MM0``.
    """
    bl = blocks(A)
    if mat24.bw24(bl[0]) != 3 or mat24.bw24(bl[1]) != 1:
        return MM0()
    blist = [x for x in range(24) if bl[0] & (1 << x)]
    dlist = sorted([ (-A[x,x], x) for x in blist])
    src = [x[1] for x in dlist]
    return MM0('p', AutPL(0, zip(src, [0,1,2]), 0)) 




def beautify_signs_size2_3(A):
    """Adjust sign of matrix with 2 times 2 and 3 times 3 blocks

    If matrix ``A`` contains 2 times 2 blocks, 3 times 3 blocks 
    and some  diagonal entries then the function returns an element 
    of the  monster that changes all off-diagonal signs to 
    reaonable values. 
    That element is returned as an instance of class ``MM0``.
    
    Otherwise the function returns the neutral element in
    class ``MM0``.
    """
    bl = blocks(A)
    if not 2 <= mat24.bw24(bl[0]) <= 3:
        return MM0()
    lplus3, lminus3, l2 = [], [], []
    for  b in bl:
        blist = [x for x in range(24) if b & (1 << x)]
        if len(blist) == 2:
            i0, i1 = blist
            l2.append((i0, i1))
        if len(blist) == 3:
            i0, i1, i2 = blist
            sign = A[i0,i1] * A[i0,i2] * A[i1,i2]
            ls = lminus3 if sign < 0 else lplus3
            l = [(i0, i1), (i0, i2), (i1, i2)]
            ls += l
    try:
        y = change_signs_A(A, lplus3 + l2, lminus3) 
        return MM0('y', y)
    except ValueError:
        try: 
            y = change_signs_A(A, lplus3 , lminus3 + l2) 
            return MM0('y', y)
        except ValueError:
            try:
                y = change_signs_A(A, l2) 
                return MM0('y', y)
            except ValueError:
                try:
                    y = change_signs_A(A, [], l2) 
                    return MM0('y', y)
                except ValueError:
                    # No improvement found
                    return MM0()

def beautify_block_size3(A):
    """Try to order a matrix with a 3 times 3 block

    If matrix ``A`` has a 3 times 3 block in rows 0, 1, and 2, then
    the function checks the diagonal entries with index >= 3. It 
    looks for a set of equal diagonal entries such that the union
    of the indices of these entries with the set [1,2,3] is an octad.
    If this is the case then the function returns an element ``g`` 
    that that moves this octad to the first 8 rows and columns.

    In case of success that element is returned as an instance of 
    class ``MM0``. Otherwise the neutral element in class ``MM0``
    is returned.    
    """
    bl = blocks(A)
    if mat24.bw24(bl[0]) != 3:
        return MM0()
    blist = [x for x in range(24) if bl[0] & (1 << x)]
    if (blist != [0,1,2]): 
        return MM0()
    d = defaultdict(list)
    for b in bl:
        if b & (b-1) == 0:
            index = b.bit_length() - 1
            d[A[index, index]].append(index)
    for lst in d.values():
        if len(lst) < 8:
            try:
                gc = GCode(blist + lst).bit_list
                assert len(gc) == 8
            except:
                continue
            isect = set(blist) & set(gc)
            #print("isect", isect)
            if len(isect) == 3:
                src0 = [x for x in gc if x in isect]
                src1 = [x for x in gc if not x in isect]
                dest = range(6)
                pi = AutPL(0, zip(src0+src1[:3], dest), 0)
                return MM0(pi)
            if len(isect) == 1:
                src = gc
                dest = list(range(8))
                src0 = [x for x in gc if not x in isect]
                src1 = [x for x in gc if x in isect]
                src2 = [x for x in blist if x not in src1]
                dest = list(range(10))
                for i in range(1000):
                    shuffle(src0)
                    shuffle(src2)
                    try:
                        pi = AutPL(0, zip(src0+src1+src2, dest))
                        return MM0(pi)
                    except:
                        pass
    return MM0()
        


#################################################################
# Dealing with one largest block of size 2
#################################################################

def beautify_block_size2(A):
    """Try to order a matrix with 2 times 2 block

    If matrix ``A`` has  2 times 2 blocks then the function checks 
    the diagonal entries. It looks for a set of equal diagonal entries 
    such that the union of the indices an octad (up to a syndrome of
    length at most 3).
    If this is the case then the function returns an element ``g`` 
    that that moves this octad to the first 8 rows and columns.

    In case of success that element is returned as an instance of 
    class ``MM0``. Otherwise the neutral element in class ``MM0``
    is returned.    
    """
    bl = blocks(A)
    if mat24.bw24(bl[0]) != 2:
        return MM0()
    blist = [x for x in range(24) if bl[0] & (1 << x)]
    d = defaultdict(list)
    for b in bl:
        if b & (b-1) == 0:
            index = b.bit_length() - 1
            d[A[index, index]].append(index)
    for lst in d.values():
        if len(lst) <= 8:
            try:
                #print(lst)
                gc = GCode(blist + lst).bit_list
                #print("GC =", gc)
            except:
                 continue
            if set(blist).issubset(gc) and len(gc) == 8:
                src0 = [x for x in gc if x in blist]
                src1 = [x for x in gc if x not in blist]
                dest = list(range(8))
                for i in range(1000):
                    shuffle(src1)
                    try:
                        pi = AutPL(0, zip(src0+src1, dest), 0)
                        return MM0(pi)
                    except:
                        pass
            if len(set(blist) & set(gc)) == 0 and len(gc) == 8:
                src0 = gc
                src1 = blist
                dest = list(range(10))
                for i in range(10000):
                    shuffle(src0)
                    try:
                        pi = AutPL(0, zip(src0+src1, dest), 0)
                        return MM0(pi)
                    except:
                        pass
    return MM0()



#################################################################
# Dealing with a diagonal matrix
#################################################################

def beautify_diagonal_matrix(A):
    """Try to order a diagonal matrix

    If matrix ``A`` is diagonal then the function looks for a set 
    of equal diagonal entries that is close to an octad or to a
    dodecad.  If this is the case then the function returns an 
    element ``g`` that that moves this Golay code word to a
    standard position.

    In case of success that element is returned as an instance of 
    class ``MM0``. Otherwise the neutral element in class ``MM0``
    is returned.    
    """
    bl = blocks(A)
    if mat24.bw24(bl[0]) != 1:
        return MM0()
    d = defaultdict(list)
    for b in bl:
        if b & (b-1) == 0:
            index = b.bit_length() - 1
            d[A[index, index]].append(index)
    singleton = None
    for lst in d.values():
         if len(lst) == 1: singleton = lst
    for lst in d.values():
        dest = list(range(8))
        if len(lst) == 8:
            #print(lst)
            for i in range(10000):
                shuffle(lst)
                try:
                    pi = AutPL(0, zip(lst, dest), 0)
                    #print("yeah")
                    return MM0(pi)
                except:
                    continue
        elif len(lst) == 11 and singleton:
            print("dodecad")
            DOD = 0xeee111
            dest = [i for i in range(24) if (1 << i) & DOD]
            src = singleton + lst
            #print(dest)
            #print(src)
            #return MM0()
            pi = mat24.perm_from_dodecads(src, dest)
            return MM0('p', pi)
        elif singleton:
            pi = pi = AutPL(0, zip(singleton, [0]), 0)
            return MM0('p', pi)            
    return MM0()




#################################################################
# Beautify signs in the special case 10B
#################################################################



def beautify_case_10B(A):
    l = sorted([(d, i) for i, d in enumerate(np.diagonal(A[:8, :8]))])
    seq = [l[j][1] for j in range(6)]
    return MM0('p', AutPL(0, zip(seq,range(6)), 0))


def beautify_signs_10B(A):
    sign_list = []
    def sgn(i0, i1, s):
        t = A[i0, i1]
        if t*s: sign_list.append((i0, i1, t*s < 0))
    for i in range(0,4): sgn(0, i+4, +1)
    for i in range(1,4): sgn(4, i, -1)
    #sgn(0, 4, 1)
    try:
        y = solve_gcode_diag(sign_list)
        print("Beautification found in case 10B")
        return MM0('y', y)
    except:
        print("Beautification failed in case 10B")
        return MM0()     



def beautify_after_signs_10B(A):
    for i in (5, 6):
        if A[1,i] < 0:
            src = [7,5,6] if i == 6 else [6,7,5]
            src = [0,1,2,3,4] + src 
            dest = list(range(7))
            #return MM0()
            return MM0('p', AutPL(0, zip(src, dest), 0))
    return MM0()     




def beautify_10B(axis):
    """Beautify an axis of type 10B"""
    axis.apply(beautify_case_10B)
    axis.apply(beautify_signs_10B)
    axis.apply(beautify_after_signs_10B)
    return axis


#################################################################
# Permuting axes
#################################################################



pi_dict = {
 "4B": zip(range(6), [2,3,0,1,4,5]),
 "6A": zip([0,1,2],[1,2,3]),
 "6C": zip([0,1,2],[1,2,3]),
 "10A": zip([0,1],[2,3]),
}

def beautify_permutation(axis):
    if axis.g_class in pi_dict:
        pi = AutPL(0, pi_dict[axis.g_class], 0)
        axis *= MM0('p', pi)


#################################################################
# Beautifying the x part of the axis
#################################################################

positive_blocks = {
   "4B": [([0,1,4,5,6,7], [0,1,4,5,6,7])],
   "6C": [([0,4,5,6,7], [0,4,5,6,7]), (0, range(4,24))],
   "10A": [(0 , range(4,24))],
   "10B": [(range(8,24), range(8,24))],
  "6A": [(0, range(4,24))],
#   "12C": [(range(7,8), range(8,16)), (range(0,4), range(4,8))],
  "4C": [(8, range(1,8))]
}


negative_blocks = {
  "8B": [(1, range(2,23))],
  "4A": [(0, range(1,23))],
   "12C": [(range(0,8), range(0,8))],
}



def map_cooctad(i):
    """Return a certain perumtation in the Mathieu group M_24

    The function returns a permutation in the  Mathieu group 
    ``M_24`` that fixes the octad (0,1,...,7) as a set, and the
    element 8, and that maps i to 9.

    The permutation is returned as an instance of class AutPL
    """
    if not 9 <= i < 24:
        raise ValueError("Bad paramter for function map_octad()")
    if i == 9:
        return AutPL(0)
    s0 = [0] + Cocode([0, 8, 9, i, i^1]).syndrome_list()
    s1 = [x for x in range(8) if not x in s0]
    src = s0 + s1[:1] + [8,i]
    dest = [0,1,2,3,4,8,9]
    return AutPL(0, zip(src,dest), 0)

def trim_cooctad(A):
    pi =  AutPL(0)  
    for i in range(9,24):
         if A[8, i]:
             pi = map_cooctad(i)
    return MM0('p', pi)

trim_cooctad_required = ["6F", "4C"]


#################################################################
#  Helper class for beautifying 2A axes
#################################################################


class Axis:
    def __init__(self, g, g_class):
        self.g = MM0(g)
        self.g_class = g_class
        self.v = v_axis * self.g
        self.v_opp = v_axis_opp * self.g
        self.g_beautify = MM0()

    def g_all(self):
        return self.g * self.g_beautify

    def __mul__(self, g):
        g = MM0(g)
        self.v *= g
        self.v_opp *= g
        self.g_beautify *= g
        return self

    @staticmethod
    def A_int(A):
        Afloat = 256 * np.copy(A)
        Aint = np.array(Afloat, dtype = np.int32)
        assert (Aint == Afloat).all(), Afloat
        return Aint

    def A(self):
        return self.A_int(self.v["A"])

    def A_opp(self):
        return self.A_int(self.v_opp["A"])

    def blocks(self):
        return blocks(self.A())

    def beautifier(self):
        self.g_beautify.reduce()
        return self.g_beautify.raw_str()
        
    def apply(self, f):
        self *= f(self.A())
        return self

    def axis_type(self, i):
        v15 = v_axis15 * self.g_all()
        return v15.axis_type(i)

    def Ax(self):
        return self.A_int((self.v * MM0('t',1))["A"])

    def reflection(self):
        return g_axis**self.g_all()        

    def trim_axis_type(self):
        t = self.axis_type(0)
        if self.axis_type(2) != t:
            self *= MM0('d', 0x800)
        if self.axis_type(1) != t: 
            Ax = self.Ax()
            if t in trim_cooctad_required:
                self *= trim_cooctad(Ax)
                Ax = self.Ax()
            if t in positive_blocks:
                d = change_signs_A(Ax, positive_blocks[t], ignore=True)
                #if (d): print("d =", d)
                self *= MM0('x', d)
            elif t in negative_blocks:
                d = change_signs_A(Ax, [], negative_blocks[t], ignore=True)
                if (d): print("d =", d)
                self *= MM0('x', d)
            else:
                self *= make_blocks_positive(Ax, tag = 'x')
        return self


#################################################################
# The beautifying function
#################################################################



def beautify(axis, verbose = 0):
    if axis.g == MM0():
        return axis
    bl = blocks(axis.A())
    if verbose:
         print("Blocksizes", [bit_count(x) for x in bl])
    bl_size = bit_count(bl[0])
    if max(bit_count(b) for b in bl) > 3:
        if bl_size > 5 and verbose: 
            print("vector type:", mat24.vect_type(bl[0]))
        axis.apply(make_blocks_positive)
        if len(bl) > 1:
            axis.apply(beautify_large_block)
        else:
            axis.apply(beautify_single_block)
        bl = axis.blocks()
        if bit_count(bl[0]) == 8 and bit_count(bl[1]) == 1:
            axis = beautify_10B(axis)
    elif bl_size == 3:
        axis.apply(beautify_large_block)
        axis.apply(sort_single_size3_block)
        axis.apply(beautify_block_size3)
        axis.apply(beautify_signs_size2_3)
    elif bl_size == 2:
        axis.apply(beautify_block_size2)
        axis.apply(beautify_signs_size2_3)
    elif bl_size == 1:
        axis.apply(beautify_diagonal_matrix)
    else:
        raise ValueError("Beautifying of matrix A failed")
    beautify_permutation(axis)
    axis.trim_axis_type()
    return axis




def compute_beautifiers(g_strings, class_names = None, verbose = 0):
    v =  MMVectorCRT(20, "I", 3, 2)
    g_data = []
    for i, g in enumerate(g_strings):
        axis = Axis(g, class_names[i] if class_names else None)
        if verbose:
            print("")
            if class_names: print("class", axis.g_class) 
            else: print("case", i)
        beautify(axis, verbose)
        #beautify_permutation(axis)
        g1, A = axis.beautifier(), axis.A()
        g_data.append(g1)
        if verbose:
            print("g =", MM0(g1))
            #print("adjacent blocks:", adjacent_blocks(A))
            print(A)
    return g_data


        



if __name__ == "__main__":
    from mmgroup.tests.test_axes.get_sample_axes import import_sample_axes
    axes = import_sample_axes()
    data = compute_beautifiers(axes.g_strings, axes.g_classes, verbose = 1)
    #print(data)

