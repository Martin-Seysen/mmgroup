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
import numpy as np
from operator import __or__
from functools import reduce

sys.path.append(r"C:\Data\projects\MonsterGit\src")

from mmgroup import MM0, MMV, MMVector, Cocode, XLeech2, Parity, PLoop
from mmgroup.mm_crt_space import MMVectorCRT 
from mmgroup import mat24, GcVector, AutPL, Parity, GCode
from mmgroup.clifford12 import bitmatrix64_solve_equation




#################################################################
# Transitive closure and identifying block of matrix A
#################################################################

def transitive_closure(m):
    """Return transitive closure of symmetric square matrix ``m``

    The function returns a 0-1 matrix ``c`` of the same size as
    matrix ``m`` if ``c[i,j] = 1`` iff there is a path of nonzero
    entries 

    ``m[k[0],k[1]], m[k[1],k[2]], ... , m[k[l-1],k[k]],``

    with ``k[0] = i, k[l] = j``. Diagoanl entrie of ``c`` are set
    to 1.
    """
    m = np.array(m != 0, dtype = np.uint32)
    assert len(m.shape) == 2 and m.shape[0] ==  m.shape[1]
    m_old = np.eye(m.shape[0], dtype = np.uint32)
    m = m | m_old
    while (m != m_old).any():
        m_old = m
        m = (m @ m) != 0
    return np.array(m, dtype = np.uint32)

def bit_count(x):
    return bin(x).count("1")

def block_sort_key(x):
    return -bit_count(x), x & -x, x 

def blocks(m):
    """Return the list of blocks of a symmetric square matrix ``m``.

    Here a block a set of indices of matrix ``m`` which are in the
    same equivalence class in the transitive closure of ``m``.

    Each block in the returned list is encoded as an integer with 
    bit ``i`` of the integer set iff index ``i`` is in that block.

    Blocks are ordered by size with larger block occuring first.
    """
    m = np.array(m != 0, dtype = np.uint32)
    m = np.array(m | m.T, dtype = np.uint32)
    m = transitive_closure(m)
    s = set()
    for row in m:
        d = int(sum(x << i for i, x in enumerate(row)))
        s.add(d)
    data = list(s)
    return sorted(data, key = block_sort_key)



def adjacent_blocks(a):
    """Return adjacent blocks of a symmetric square matrix ``m``.

    If matrix ``m`` may be split into ``n`` adjacent blocks then  
    the function returns a list ``i`` of lenggth ``n+1`` such that
    all nonzero entries of ``m`` occur in one of the submatrices

    ``m[i[0]:i[1], i[0]:i[1]], ...  m[i[n]:i[n+1], i[n]:i[n+1]]``

    Then ``i[0] = 0`` and `i[n+1]`` is equal to the size of the
    square matrix ``m``.  
    """
    bl = blocks(a)
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


def format_eigen_values(a):
   eigen =  purge_diag(np.linalg.eigvals(a))
   eigenv = sorted(eigen.items(), key = lambda x: x[0])
   data = [("%d^%d" if isinstance(k, int) else "%.2f^%d") % 
       (k,v) for k, v in eigenv ]
   return "(" + ", ".join(data) + ")"
 



def block_eigenvalues(a):
    ad = adjacent_blocks(a)
    if len(ad) > 24:
        return format_eigen_values(a)
    eigen = []
    exp = []
    for i in range(len(ad) - 1):
        m = a[ad[i]:ad[i+1], ad[i]:ad[i+1]] 
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


def beautify_large_block(bl):
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



def make_blocks_positive(A, first_row_only = False):
    bl = blocks(A)   
    bl_list = [x for x in range(24) if (1 << x) & bl[0]] 
    equations = []
    for i in range(24):
        for j in range(24):
            if i == j or A[i,j] == 0:
                continue
            if i in bl_list and j in bl_list:
                s = A[i,j] < 0
                if not first_row_only or i == bl_list[0]:
                    equations.append((i, j,  s))
            else:
                equations.append((i, j,  0))
    try:
        y = solve_gcode_diag(equations)             
        return MM0('y', y)
    except:
        if not first_row_only:
            try:
                return make_blocks_positive(A, True)
            except:
                pass
        return MM0()



#################################################################
# Dealing with one largest block of size 3
#################################################################



def beautify_signs_size3(A):
    bl = blocks(A)
    if mat24.bw24(bl[0]) != 3:
        return MM0()
    sign = A[0,1] * A[0,2] * A[1,2]
    if sign > 0:
        l = [(i, j, A[i,j] < 0)  for i,j in [(0,1),(0,2),(0,3)]]
        y = solve_gcode_diag(l)
        if y:
            return MM0('y', y)
    return MM0()



def beautify_signs_size2_3(A):
    bl = blocks(A)
    if not 2 <= mat24.bw24(bl[0]) <= 3:
        return MM0()
    data2, data3 = [], []
    for  b in bl:
        blist = [x for x in range(24) if b & (1 << x)]
        if len(blist) == 2:
            i0, i1 = blist
            data2.append((i0, i1, A[i0,i1] < 0))
        if len(blist) == 3:
            i0, i1, i2 = blist
            sign = A[i0,i1] * A[i0,i2] * A[i1,i2]
            for j0, j1 in [(i0,i1), (i0,i2), (i1,i2)]:
                s = sign * A[j0,j1] < 0
                data3.append((j0, j1, s < 0))
    try:
        y = solve_gcode_diag(data2 + data3)
        return MM0('y', y)
    except:
        try:
            data2n = [(i0,i1,1-x) for i0,i1,x in data2]
            y = solve_gcode_diag(data2n + data3)
            return MM0('y', y)
        except:
            raise ValueError("Solve")
            pass
    return MM0()

def beautify_block_size3(A):
    bl = blocks(A)
    if mat24.bw24(bl[0]) != 3:
        return MM0()
    blist = [x for x in range(24) if bl[0] & (1 << x)]
    assert blist == [0,1,2], (blist, A)
    d = defaultdict(list)
    for b in bl:
        if b & (b-1) == 0:
            index = b.bit_length() - 1
            d[A[index, index]].append(index)
    for lst in d.values():
        if len(lst) < 8:
            try:
                #print(lst)
                gc = GCode(blist + lst)
                #print("GC =", gc.bit_list)
            except:
                #print("WTF")
                raise
            len_isect = len(set(blist) & set(gc.bit_list))
            #print("isect", len_isect)
            if len_isect == 3:
                src = gc.bit_list
                dest = list(range(8))
                return MM0(AutPL(0, zip(src, dest), 0))
            if len_isect == 1:
                src0 = gc.bit_list[1:]
                src1 = gc.bit_list[:1]
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
    bl = blocks(A)
    if mat24.bw24(bl[0]) != 2:
        return MM0()
    blist = [x for x in range(24) if bl[0] & (1 << x)]
    #assert blist == [0,1], blist
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
                 return MM0()
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
            print(dest)
            print(src)
            #return MM0()
            pi = mat24.perm_from_dodecads(src, dest)
            return MM0('p', pi)
        elif singleton:
            pi = pi = AutPL(0, zip(singleton, [0]), 0)
            return MM0('p', pi)            
    return MM0()




#################################################################
# Beautify signs in the secial case 10N
#################################################################

def find_five_cycle(A):
   """Try to find a cycle of positive off-diagonal elemnts in submatrix

   Let ``M `` be the subatrix ``A[1:6, 1:6]`` of the symmetric matrix 
   ``A``. The function entries to find a cycle of length 5 such that
   ``A[i0, i1] > 0`` for any adjacent entries ``i0, i1`` of the cycle.

   If such a cycle exists then the function returns that cycle as a list
   of integers between 1 and 6, correspondig to the indices of ``A``.
   If no such cycle exist the the function returns the empty list.
 
   The function checks the following condition (and returns ``[]`` if
   any of them is not satisfied):

   * ``A[:6, 6:]`` must be zero.

   * The diagonal of ``A[:6]`` must be nonnegative.

   * A[0,:] must be nonnegative.

   * ``A[i, 1:6], 1 <= i < 6`` must have precisely two (off-diagonal)
     negative entries.   
   """
   if min(A[0]) < 0: 
       return []
   for i in range(1,6):
       x = sum(x < 0 for x in A[i])
       if (x != 2) or A[i,i] < 0:
           return []
   if max(abs(A[:6, 6:].ravel())) > 0:
       return []
   l = [1]
   for i in range(5):
       for j in range(1,6):
           if A[j, l[-1]] > 0 and j not in l:
               l.append(j)
               break
   if len(l) < 5 or A[l[0],l[-1]] < 0:
       return []
   return l



def beautify_case_10B(A):
    bl = blocks(A)
    if bl[:3] != [63, 192, 256]:
        return MM0()
    c = find_five_cycle(A)
    if len(c) != 5:
        return MM0()
    for j in  [0, 1]:
        pi_src = [0] + c + [6+j, 7-j]
        pi_dest = list(range(8))
        try:
           pi = AutPL(0, zip(pi_src, pi_dest), False)
        except ValueError:
           continue
        return MM0('p', pi)
    return MM0()     


#################################################################
# Enumarating the 2A axes
#################################################################



def get_A(w):
    Afloat = 128 * w["A"]
    A = np.array(Afloat, dtype = np.int32)
    assert (A == Afloat).all()
    return A



def beautify(v, verbose = 0):
    g = MM0()
    A = get_A(v)
    bl = blocks(A)
    if verbose:
         print("Blocksizes", [bit_count(x) for x in bl])
    bl_size = bit_count(bl[0])
    if bl_size > 3:
        if bl_size > 5 and verbose: 
            print("vector type:", mat24.vect_type(bl[0]))
        g1 = beautify_large_block(bl[0])
        v *= g1
        g *= g1
        A = get_A(v)
        g2 = make_blocks_positive(A)
        v *= g2
        g *= g2
        A = get_A(v)
        g3 = beautify_case_10B(A)
        v *= g3
        g *= g3
        return g.reduce(), get_A(v)
    elif bl_size == 3:
        g0 = beautify_large_block(bl[0])
        v *= g0
        g *= g0
        A = get_A(v)
        g1 = beautify_signs_size3(A)
        v *= g1
        g *= g1
        A = get_A(v)        
        g2 = beautify_block_size3(A)
        v *= g2
        g *= g2
        A = get_A(v)
        g3 = beautify_signs_size2_3(A)
        v *= g3
        g *= g3
        return g.reduce(), get_A(v)
    elif bl_size == 2:
        g1 = beautify_block_size2(A)
        v *= g1
        A = get_A(v)
        g2 = beautify_signs_size2_3(A)
        v *= g2
        return (g1 * g2).reduce(), get_A(v)
    elif bl_size == 1:
        g1 = beautify_diagonal_matrix(A)
        v *= g1
        return (g1).reduce(), get_A(v)
    else:
        raise ValueError("Beautifying of matrix A failed")



def compute_beautifiers(g_strings, class_names = None, verbose = 0):
    v =  MMVectorCRT(20, "I", 3, 2)
    g_data = []
    for i, g in enumerate(g_strings):
        if verbose:
            print("")
            if class_names: print("class", class_names[i]) 
            else: print("case", i)
        if  MM0(g) != MM0():
            g1, A = beautify(v * MM0(g), verbose)
        else:
            g1, A =  MM0(), get_A(v)
        g_data.append(g1.raw_str())
        if verbose:
            print("g =", g1)
            print("adjacent blocks:", adjacent_blocks(A))
            print(A)
    return g_data


if __name__ == "__main__":
    from mmgroup.tests.test_axes.get_sample_axes import import_sample_axes
    axes = import_sample_axes()
    data = compute_beautifiers(axes.g_strings, axes.g_classes, verbose = 1)
    #print(data)

