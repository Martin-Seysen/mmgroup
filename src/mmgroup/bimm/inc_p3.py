r"""The module implements the projective plane over :math:`\mathbb{F}_3`


Module ``inc_p3`` implements the projective plane ``P3`` over
:math:`\mathbb{F}_3` and its automorphism group. Points and lines
in ``P3`` are implemented as instances of class ``P3_node``,
automorphisms of ``P3`` are implemented as instances of 
class ``AutP3``. 

We mainly work with the incidence graph containing the points
and the lines of ``P3`` as vertices.
The union of the set of the points and of the set of the lines 
in ``P3`` is called the set of the *nodes* of the projective 
plane ``P3``.
"""
import os
import sys
from numbers import Integral
from random import randint, sample, choice
from functools import reduce
from operator import __and__, __or__

import numpy as np


import mmgroup
from mmgroup.clifford12 import uint64_to_bitlist
from mmgroup.clifford12 import uint64_low_bit
from mmgroup.clifford12 import uint64_to_bitarray
from mmgroup.clifford12 import uint64_bit_weight
from mmgroup.structures.abstract_group import AbstractGroup
from mmgroup.structures.abstract_group import AbstractGroupWord
from mmgroup.structures.abstract_group import singleton



# The following code block is just for generating documentation
# with readthedocs and can be ignored.
try:
    # A stupid way to circumvent the mockup process for readthedocs
    assert uint64_to_bitlist(3) == [0,1] # check if function is correct
    #assert False
except:
    # Reprogram (slow) low-level functions if fast versions not present
    def uint64_to_bitlist(n):
        n = int(n)
        return [i for i in range(64) if (n >> i) & 1]
    def uint64_bit_weight(n):
        n = int(n)
        return sum([(n >> i) & 1 for i in range(64)])
    def uint64_low_bit(n):
        n = int(n) & 0xffffffffffffffff
        if n == 0:
            return 64
        return (n & -n).bit_length() - 1

# Some error messages
ERR_PROJ = "Mapping does not preserve the projective plane P3"
ERR_UNIQUE = "Mapping is underdetermined in the projective plane P3"
ERR_DUPL = "Duplicate entry in mapping of projective plane P3"
ERR_PL_ALL = "P3 nodes in %s must all be points or all lines"
ERR_P_ALL = "P3 nodes in %s must all be points"


#####################################################################
# Names of nodes in the projective plane P3
#####################################################################

# Dictionary ``P3_OBJ`` maps integers and string corrsponding to
# number or names of nodes in the projective
P3_OBJ = dict(zip(range(26), range(26)))
for x in range(26):
    P3_OBJ[str(x)] = x
for x in range(13):
    P3_OBJ['P'+str(x)] = x
    P3_OBJ['L'+str(x)] = x + 13


    
def p3_obj(obj):
    r"""Convert python object ``obj`` to the number of a P3 node"""
    if isinstance(obj, P3_node):
       return obj.ord
    try:
       return P3_OBJ[obj]
    except KeyError:
       if isinstance(obj, Integral):
           raise IndexError("Number of P3 node out of range")
       elif isinstance(obj, str):
           s = " '%s'" % obj if len(obj) < 11 else ""
           err = "Cannot convert string%s to P3 node"
           raise ValueError(err % s)
       else:
           err = "Cannot convert %s object to P3 node"
           raise TypeError(err % str(type(obj)))


def p3_list(obj):
    r"""Convert python object ``obj`` to list of P3 nodes"""
    if isinstance(obj, str):
        s = [x.strip() for x in obj.split(',') if not x.isspace()]
        return [p3_obj(x) for x in s]
    else:
        return list(map(p3_obj, obj))
 




class P3_node:
    r"""Models a point or a line the projective plane P3.

    We number the 13 points in the projective plane ``P3`` over
    :math:`\mathbb{F}_3` from 0 to 12, and the 13 lines from
    13 to 25. Then a point with number ``i`` and a line with
    number ``j`` are incident if
    :math:`i + j \equiv 0, 1, 3,` or :math:`9 \pmod{13}`.

    Some strings are also accepted as a description of a point
    or a line in ``P3``. The 13 points may be denoted as
    'P0',...,'P12', and the the 13 lines may be denoted as
    'L0',...,'L12'.

    The names 'a', 'b1', 'b2', 'b3', 'c1', 'c2', 'c3', etc.
    refer to the embedding of the :math:`Y_{555}` graph into
    the projective plane ``P3``, as described in the documentation
    of the application. For background, see also :cite:`CNS88`,
    :cite:`Far12`.

    :param obj: An integer, a string, or an instance of
                class ``AutP3`` describing a point or a line 
                in the projective plane ``P3``.
    """
    slots = ['_ord']
    def __init__(self, obj):
        self._ord = int(p3_obj(obj)) 
    def __str__(self):
        t = "point" if self._ord < 13 else "line"
        return "P3<%s %d>" % (t, self._ord % 13)
    def __eq__(self, other):
        return isinstance(other, P3_node) and self._ord == other._ord
    def __ne__(self, other):
        return not self.__eq__(other)
    def __mul__(self, other):
        if isinstance(other, AutP3):
            o = self._ord
            if o < 13:
                return P3_node(other.perm[o])
            else:
                p1, p2 = INCIDENCE_LISTS[o, :2]
                im1, im2 = other.perm[p1], other.perm[p2] 
                return P3_incidence(im1, im2) 
        else:
            err = "Cannot multiply class P3_node object with %s object"
            raise ValueError(err % type(other))
    @property
    def ord(self):
        r"""Return internal number of instance of class ``P3_node``"""
        return self._ord
    def name(self):
        r"""Return the name of the ``P3`` node in standard notation"""
        q, r = divmod(self._ord, 13)
        return "PL"[q] + str(r)
    def y_name(self):
        r"""Return the name of the ``P3`` node in Y_555 notation"""
        return Y_NAMES[self._ord]

#####################################################################
# Elementary geometry in projective plane P3
#####################################################################

# Bit map containing all 26 nodes of the projective plane P3
ALL_BITS = 0x3ffffff

# INCIDENCES[x] is the list of integers. Bit y of INCIDENCES[x] 
# is set if node y is incident with node x.
INCIDENCES = np.zeros(26, dtype = np.uint32)
# INCIDENCES_LISTS[x] is the array of the 4 nodes incident with x 
INCIDENCE_LISTS = np.zeros((26,4), dtype = np.uint32)
for x in range(13):
    blist = sum(1 << ((p - x) % 13) for p in (0, 1, 3, 9))
    INCIDENCES[x] = blist << 13
    INCIDENCES[x + 13] = blist
    INCIDENCE_LISTS[x] = uint64_to_bitlist(INCIDENCES[x])
    INCIDENCE_LISTS[x + 13] = uint64_to_bitlist(INCIDENCES[x + 13])



def P3_incidences(*x):
    r"""Return list of P3 nodes incident with given P3 nodes

    Here each argument of the function is a list of nodes of P3
    describing a set :math:`S_i` of nodes (i.e. points or lines) of
    P3. An entry of such a list may be anything that is accepted by
    the constructor of class ``P3_node``. An integer argument is
    interpreted as a singleton, i.e. a set of size 1. A comma-
    separated string of names of P3 nodes is accepted as a set
    of P3 nodes.

    The function returns the sorted list of P3 nodes (i.e instances of
    class ``AutP3`` that are incident with at least one node in each
    set :math:`S_i` and not contained in any of the sets :math:`S_i`.
    """
    nodes = ALL_BITS
    no_nodes = 0
    for l in x:
        if isinstance(l, Integral):
            nodes &= INCIDENCES[l]
            no_nodes |= 1 << l
        else:
            l1 = p3_list(l)
            nodes &= reduce(__or__, [INCIDENCES[p] for p in l1], 0)
            no_nodes |= reduce(__or__, [1 << p for p in l1], 0)
    nodes = uint64_to_bitlist(int(nodes) & ~no_nodes)
    return [P3_node(x) for x in nodes]


def P3_incidence(*x): 
    r"""Return (unique) P3 node incident with given P3 nodes

    Here each argument describes a P3 node (i.e. a point
    or a line); it may be anything accepted by the constructor
    of class ``P3_node``.

    If there is a unique P3 node incident with all these P3
    nodes then the function returns that node as an instance of
    class ``P3_node``.

    Otherwise the function raises ValueError.

    This function is a simplified version of function
    ``P3_incidences``. Its typical use case is to find a line
    through two points or the intersection point of two lines.
    """ 
    a = reduce(__and__, [INCIDENCES[p3_obj(p)] for p in x], ALL_BITS)
    if uint64_bit_weight(a) == 1:
        return  P3_node(uint64_low_bit(a))
    if (a):
        s = "Incident node in function incidence() is not unique"
    else:
        s = "No incident node found in function incidence()"
    raise ValueError(s)
    


    

def _remaining_nodes(x1, x2):
    r"""Internal version of function` `P3_remaining_nodes``

    Input and operation is as in function ``P3_remaining_nodes``.
    But the result is returned as a list of integers.
    """
    blist = uint64_to_bitlist(
        int(INCIDENCES[p3_obj(x1)]) & int(INCIDENCES[p3_obj(x2)])
    )
    if len(blist) == 1:
        rem = int(INCIDENCES[blist[0]]) & ~((1 << x1) | (1 << x2))
        return uint64_to_bitlist(rem)
    if len(blist):
        s = "Arguments in P3_remaining_nodes() must be differnt"
    else:
        s = ERR_PL_ALL % 'P3_remaining_nodes()'
    raise ValueError(s)

def P3_remaining_nodes(x1, x2):
    r"""Complete points on a line or lines intersecting in a point

    If arguments ``x1, x2`` are different points or lines in P3
    then the function returns the list of the two remaining points
    on the line through ``x1`` and ``x2``, or the list of the two
    remaining lines containing the intersection point of ``x1`` and
    ``x2``, respectively. Otherwise the function raises ValueError.

    Arguments ``x1, x2`` may by anything accepted by the constructor
    of class ``P3_node``. The result is returned as list of two     
    instances of class ``P3_node``.
    """
    return [P3_node(x) for x in remaining_nodes(x1, x2)]



def find_cross(points):
    r"""Find quadruple of non-collinear points in list of points

    Let ``points`` be a list of different points in the projective
    plane P3. Define a **cross** to be a set of four points in P3,
    with no three points in that set being collinear.

    The function returns a cross in ``points`` as a list of
    for points in ``points``. It returns ``None`` if the set
    ``points`` does not contain any cross. Duplicate entries
    in the list ``points`` are illegal.

    Any set of 6 or more points contains a cross.
    """
    def common(x1, x2):
        """Return list of all points on line through x1 and x2"""
        blist = uint64_to_bitlist(INCIDENCES[x1] & INCIDENCES[x2])
        assert len(blist) == 1
        return int(INCIDENCES[blist[0]])

    # Four different points x1, x2, x3, x4 form a cross if the inter-
    # section of the line through x1 and x2, and the line through
    # x3 and x4, contains a point differnt from x1, x2, x3, x4.
    if len(set(points)) < len(points):
        raise ValueError(ERR_DUPL)
    if not 0 <= min(points) <= max(points) < 13:
         raise ValueError(ERR_P_ALL % 'list')
    points = [int(x) for x in points[:6]]
    n = len(points)
    if n < 4:
        return None
    for i1, x1 in enumerate(points):
        for i2 in range(i1 + 1, n):
            x2 = points[i2]
            s12 = common(x1, x2)
            for i3 in range(i2 + 1, n):
                x3 = points[i3]
                if s12 & (1 << x3):
                    continue
                s123 = s12 | common(x1, x3) | common(x2, x3)
                for i4 in range(i3 + 1, n):
                    x4 = points[i4]
                    if (s123 & (1 << x4)) == 0:
                       return [x1, x2, x3, x4]
    return None            
 
# Store a fixed cross in variable FST_CROSS
FST_CROSS = find_cross(range(13))
assert FST_CROSS == [0, 1, 2, 5]


def find_collinear_points(points):
    r"""Find triple of collinear points in list of points

    Let ``points`` be a list of different points in the projective
    plane P3. The function tries to find a triple ``(x1, x2, x3)``
    of collinear points in the set ``points``. If such a triple
    exists then the function returns a list ``[x1, x2, x3, x4]``,
    where ``(x1, x2, x3)`` is the triple found, and ``x4`` is the
    fourth point on the line given by that triple. Here ``x4`` may
    or may not be in the set ``points``. Duplicate entries
    in the list ``points`` are illegal.


    The function returns ``None`` if the set``points`` does not 
    contain three collinear points.

    Any set of 5 or more points contains three collinear points.
    """
    if len(set(points)) < len(points):
        raise ValueError(ERR_DUPL)
    if len(points) and not 0 <= min(points) <= max(points) < 13:
         raise ValueError(ERR_PL_ALL % 'points')
    if len(points) < 3:
        return None
    points = points[:5]
    for i1, x1 in enumerate(points):
        for x2 in points[i1 + 1:]:
            x3, x4 = _remaining_nodes(x1, x2)
            if x3 in points:
                return [int(x1), int(x2), int(x3), int(x4)]
            if x4 in points:
                return [int(x1), int(x2), int(x4), int(x3)]
    return None            
       


def P3_is_collinear(l):
    r"""Check if list of P3 nodes contains 3 collinear nodes

    Argument ``l`` of the function is a list of nodes of P3
    describing a set of nodes (i.e. of points or lines) of P3. 
    An entry of such a list may be anything that is accepted by
    the constructor of class ``P3_node``.  A comma-separated
    string of names of P3 nodes is accepted as a set of nodes.

    The function returns ``True`` if that set of nodes contains 3 
    collinear points or 3 collinear lines, and ``False`` otherwise.
    """
    bitmap = reduce(__or__, [1 << x for x in p3_list(l)], 0)
    if find_collinear_points(uint64_to_bitlist(bitmap & 0x1fff)):
        return True
    if find_collinear_points(uint64_to_bitlist(bitmap >> 13)):
        return True
    return False



def complete_cross_random(points):
    r"""Find or complete quadruple of non-collinear points

    Let ``points`` be a list of different points in the projective
    plane P3. Define a **cross** to be a set of four points in P3,
    with no three points in that set being collinear.

    The function returns a cross in  the list ``points`` as a list
    of  for points if found. If the list ``points`` does not
    contain any cross, then the function (internally) completes that 
    list with random points until a cross is found , and returns
    that cross.

    Duplicate entries in the list ``points`` are illegal.
    """
    point = [int(x) for x in points]
    if len(points) == 0:
        points = [randint(0,12)]
    cross = find_cross(points) 
    if cross:
        return cross
    points = points[:6]
    # Here ``points`` contains at most one set of 3 or 4 
    # collinear points
    line = find_collinear_points(points)
    if line:
        others = list(set(points) - set(line))
        if len(others) == 0:
            y1 = choice(tuple(set(range(13)) - set(line)))
        else:
            assert len(others) == 1
            y1 = others.pop()
        y2 = choice(_remaining_nodes(y1, line[2]))
        return [line[0], line[1], y1, y2]
    # Here ``point`` contains at most 3 (non-collinear) points
    # Fill ``points`` with random points up to length at least 2
    others = set(range(13)) - set(points)
    if len(points) < 2: 
        points += sample(tuple(others), 2 - len(points))
        others = others - set(points)
    # Add a 3rd non-collinear (random) point if not present
    others = others - set(_remaining_nodes(*points[:2]))
    if len(points) < 3:
        points.append(choice(tuple(others)))
        others.remove(points[2])
    # Add a 4-th non-collinear (random) point
    others = others - set(_remaining_nodes(points[0], points[2]))
    others = others - set(_remaining_nodes(points[1], points[2]))
    points.append(choice(tuple(others)))
    return points





def cross_intersection(x11, x12, x21, x22):
    r"""Intersection of two lines, each line given by to points

    We assume that no three of the points ``x11, x12, x21, x22``
    in the projective plane P3 are collinear; otherwise the
    function raises ValueError.

    Let ``l1`` be the line through points ``x11, x12``, and let
    ``l2`` be the line through points ``x21, x22``. Let ``y``
    be the point at the intersection of ``l1`` and ``l2``.
    Let ``y1`` and ``y2`` be the remaining points of the
    four points on lines ``l1`` and ``l2``, respectively.

    The function returns the list ``[y, y1, y2]``
    """
    def remain(x1, x2):
        r"""Return bitnap of remaining points on line through x1 and x2"""
        blist = uint64_to_bitlist(INCIDENCES[x1] & INCIDENCES[x2])
        assert len(blist) == 1
        return int(INCIDENCES[blist[0]]) &~ ((1 << x1) | (1 << x2))

    s1 = remain(x11, x12)
    s2 = remain(x21, x22)
    #print("Cross intersection", x11, x12, x21, x22, hex(s1), hex(s2))
    if s1 == s2 or (s1 & s2) == 0:
        s = "Collinear points found in function cross_intersection()"
        raise ValueError(s)
    return [uint64_low_bit(s1 & s2), uint64_low_bit(s1 & ~s2), 
              uint64_low_bit(~s1 & s2)]




def map_cross(cross1, cross2):
    r"""Return unique mapping from one cross to another cross

    Let a *cross* be an (ordered) list of four non-collinear points 
    in P3 as in function ``find_cross``. If ``cross1`` and ``cross2``
    are two such crosses then there is a unique automorphism of
    P3 that maps ``cross1`` to ``cross2``. 

    The function computes that automorphism as a permutation of the
    13 points of P3. It returns a list, where entry ``i`` of that
    list is the image of point ``i``.
    """
    perm = [-1] * 13
    assert len(cross1) == len(cross2) == 4
    c1 = [x % 13 for x in cross1]
    c2 = [x % 13 for x in cross2]
    for i in range(3):
        c1 += cross_intersection(*c1[:4])
        c2 += cross_intersection(*c2[:4])
        c1[0:3] = c1[1], c1[2], c1[0]
        c2[0:3] = c2[1], c2[2], c2[0]
    for i in range(13):
        perm[c1[i]] = c2[i]
    assert min(perm) == 0 and max(perm) == 12
    return perm



def line_map_from_map(perm):
    r"""Convert mapping of points of P3 to a mapping of lines

    Let ``perm`` be a list of length 13 such that the mapping
    ``i -> perm[i]`` of points is an automorphism of the projective 
    plane P3.

    The function returns a list of length 13 containing the same
    mapping as a permutation of lines. Entry ``i`` of the returned 
    list is the image of line  ``i``. Here all line numbers are 
    reduced modulo 13.

    Due to the symmetry between the point numering and the line 
    numbering, this function may also be used to convert mapping of 
    lines of P3 to a mapping of points.
    """
    line_perm = []
    for x in range(13,26):
        p1, p2 = INCIDENCE_LISTS[x,:2]
        img = INCIDENCES[perm[p1]] & INCIDENCES[perm[p2]]
        line_perm.append(uint64_low_bit(int(img) >> 13))
    return line_perm


def map_P3_to_perm(obj1, obj2, unique = True):
    r"""Convert mapping of P3 nodes to permutation of points

    Arguments ``obj1`` and ``obj2`` are lists of integers of
    the same length that define a (partial) mapping
    ``obj1[i] -> obj2[i]`` of P3 nodes. The function tries
    to find an automorphism of P3 compatible with that mapping.
    In the current version the entries of the lists ``obj1`` 
    and ``obj2`` must all be points or all lines.

    In case of success the function returns that automorphism as 
    a permutation of the 13 points of P3. It returns a list, where 
    entry ``i`` of that list is the image of point ``i``.
    The function raises ValueError in case of failure.

    If parameter ``unique`` is ``True`` (default) then the function
    fails if the mapping ``obj1[i] -> obj2[i]`` does not extend
    to a unique automorphism of P3. Otherwise the function
    returns a (uniform distributed) random feasible automorphism. 
    """
    all_obj_by13 = [x // 13 for x in obj1 + obj2]
    line = 0
    if len(all_obj_by13):
        line = min(all_obj_by13)
        if not 0 <= line == max(all_obj_by13)  <= 1:
            raise ValueError(ERR_PL_ALL % 'mapping')
    if len(obj1) != len(obj2):
        err = "Preimage and image of mapping must have same length"
        raise ValueError(err)
    obj1 = [x % 13 for x in obj1]
    obj2 = [x % 13 for x in obj2]
    if unique:
        cross1 = find_cross(obj1) 
        if not cross1:
            raise ValueError(ERR_UNIQUE)
        #cross2 = [obj2[obj1.index(x)] for x in cross1]
        cross2 = find_cross(obj2) 
    else:
        cross1 = complete_cross_random(obj1)
        cross2 = complete_cross_random(obj2)
    perm = map_cross(cross1, cross2)
    for i, p1 in enumerate(obj1):
        if perm[p1] != obj2[i]:
            raise ValueError(ERR_PROJ)
    if line:
        perm = line_map_from_map(perm)
    return perm



def check_perm_P3(perm):
    r"""Check if a mapping of points in P3 is an automorphism

    Let a mapping of points in P3 be given by ``i-> perm[i]``,
    for ``0 <= i < 13´`. The function checks if that mapping is
    an automorphism of ``P3``. It returns ``perm`` if this is 
    the case. Otherwise it raises ValueError.
    """
    if len(perm) != 13:
        s = "Point permutation list must have length 13"
        raise ValueError(s)
    perm =  [x % 13 for x in perm]
    img_cross = [perm[i] for i in FST_CROSS]
    img_perm = map_P3_to_perm(FST_CROSS, img_cross)
    if perm != img_perm:
        raise ValueError(ERR_PROJ)
    return perm


def invert_perm_P3(perm):
    r"""Invert a permutation of points in P3

    Let a permutation of points in P3 be given by ``i-> perm[i]``,
    for ``0 <= i < 13´`. The function returns the inverse of
    that permutation.

    The resultis undefined if the permutation is not an automorphsim
    of P2.
    """
    i_perm = [None] * 13
    for i, x in enumerate(perm):
        i_perm[x] = i
    check_perm_P3(perm)
    return i_perm


def mul_perm_P3(perm1, perm2):
    r"""Multiply two permutations of points in P3

    Let a ``perm1, perm2`` be permutations of points in P3 be given 
    by ``i-> perm1[i], i-> perm2[i]``,   for ``0 <= i < 13´`. The 
    function returns the product ``perm1 * perm2`` of these two
    permutations.

    The resultis undefined if the permutation is not an automorphsim
    of P2.
    """
    perm3 = [None] * 13
    for i, x in enumerate(perm1):
        perm3[i] = int(perm2[x])
    check_perm_P3(perm3)
    return perm3



         


def P3_point_set_type(s):
    r"""Return a certain invariant of a set of points of P3

    Given a list or a set ``s`` of points of P3, the function
    returns a certain tuple depending on ``s`` that is 
    invariant under the automorphism group of P3.
    """
    assert 0 <= min(s)
    assert max(s) < 13
    bl = reduce(__or__, [1 << p for p in s], 0)
    l = [0]*5
    for i in range(13, 26):
        l[uint64_bit_weight(int(INCIDENCES[i]) & bl)] += 1
    return uint64_bit_weight(bl), tuple(l)
        

#####################################################################
# Add Y_555 names to dict P3_OBJ
#####################################################################


# Enter names of the nodes in Y_555 notation into dictionary P3_OBJ
# We start with assigning numbers to the nodes 'a', 'c1', 'c2', 'c3'.
P3_OBJ.update({'a':0, 'c1':1, 'c2':2, 'c3':11})
assert find_cross([P3_OBJ[x] for x in ('a', 'c1', 'c2', 'c3')])
def _join(args):
    r"""Add names of the nodes in Y_555 notation into dictionary P3_OBJ

    Here 'args' is a comma-sparated string of triples of alphanumeric
    names of nodes. Those triples are separated by whitespace. 

    A triple ``(x, y, z)`` means that ``z`` is the unique node
    incident with the two nodes ``x`` and ``y``. Here the names ``x``
    and ``y`` must already be keys in the dictionary ``P3_OBJ``.
    """
    triples = [x.strip() for x in args.split(',') if not x.isspace()]
    for s in triples:
        a = s.split(' ')
        if len(a) == 3: 
            P3_OBJ[a[0]] = P3_incidence(a[1], a[2])._ord 
# Enter all nodes (in Y_555 notation) into the dictionary P3_OBJ.
# Starting with the above assignment of nodes the construction is forced.
# Names of the nodes and their incidences are taken from the ATLAS.
_join('b1 a c1, b2 a c2, b3 a c3, z1 c2 c3, z2 c1 c3, z3 c1 c2')
_join('a1 b1 z1, a2 b2 z2, a3 b3 z3, c1 z2 z3, c2 z1 z3, c3 z1 z2')
_join('f1 a2 a3, f2 a1 a3, f3 a1 a2, g1 b1 f1, g2 b2 f2, g3 b3 f3')
_join('d1 g2 g3, d2 g1 g3, d3 g1 g2, e1 z1 d1, e2 z2 d2, e3 z3 d3')
_join('f a e1')



# Dictionary Y_names maps the numbers of the nodes to their Y_555 names
Y_NAMES = {}
for k, v in P3_OBJ.items():
    if isinstance(k, str) and k[0] in "abcdefgz":
         Y_NAMES[v] = k 
assert len(Y_NAMES) == len(set(Y_NAMES.values())) == 26

             



#####################################################################
# class AutP3 
#####################################################################





def p3_mapping(src = None, random = False):
    """Implements the constructor of class AutP3, standard cases"""
    if src is None:
        if not random:
            return list(range(13))
        else:
            src = {0:randint(0,12)}
    if isinstance(src, str):
        try:
            s = [x.strip() for x in src.split(',') if not x.isspace()]
            src = dict([[y.strip() for y in x.split(':')] for x in s])
        except:
            err = "Cannot evaluate string to a mapping of P3 nodes"
            raise ValueError(err)
    elif isinstance(src, zip):
        src = dict(src)
    if isinstance(src, dict):
        h1, h2 = zip(*src.items())
        return map_P3_to_perm(p3_list(h1), p3_list(h2), not random)
    s = "Cannot construct automorphism of P3 from '%s' object"     
    raise TypeError(s % type(src))




NEUTRAL_PERM_P3 = list(range(13))




class AutP3(AbstractGroupWord):
    r"""Models an automorphism of the projective plane ``P3``.

    This class models the automorphism group ``AutP3`` of the 
    projective plane ``P3`` over the field :math:`\mathbb{F}_3`.

    Elements of ``AutP3`` should be given as (partial) mappings of 
    points or lines. The standard way to describe an automorphism  
    in ``AutP3`` is a dictionary containing a partial mapping 
    of points or lines. Here the keys and the values of the
    dictionary must either all be points or all lines; they must be
    objects describing points or lines as in the constructor of
    class ``AutP3``. A mapping between points or lines is accepted
    if it extends to a unique mapping of the projective plane ``P3``. 


    :param mapping:

      Describes a mapping of points or lines in the projective
      plane ``P3`` as indicated in the table below.

    :param data:

      Additional data (optional) that describe a mapping of points 
      or lines in some special cases.
   

    .. table:: Legal types for parameter ``mapping`` in the constructor
      :widths: 25 75

      ================ ===================================================
      type             Evaluates to
      ================ ===================================================
      ``None``         Creates the neutral element (default).
  

      class ``AutP3``  A deep copy of the given automorphism in class
                       ``AutP3`` is returned. 
  
      ``dict``         Dictionary containing a mapping between points or
                       lines as described above.

      ``zip`` object   ``zip(x,y)`` is equivalent to ``dict(zip(x,y))``

      ``string`` 'r'   Then we construct a random automorphism (depending
                       on parameter ``data``) as described below. 

      ``string`` 'p'   Then ``data`` must be a list of 13 integers (taken
                       modulo 13), that describes a mapping of the 
                       13 points. 

      ``string`` 'l'   Then ``data`` must be a list of 13 integers (taken
                       modulo 13), that describes a mapping of the 
                       13 lines. 
      ================ ===================================================



    Remarks:

    If parameter ``mapping`` is the string ``'r'``, then an optional
    parameter ``data`` of type ``dict`` or ``zip`` that describes a 
    partial mapping of points or lines may follow. In this case we 
    construct a random  automorphism of ``P3`` satifying the 
    constraints of the mapping given by parameter ``data``, if present.
    Such a random automorphism  is chosen from a uniform distribution 
    of all possible cases.

    For instances ``g1`` and ``g2`` of this class,
    ``g1 * g2``  means group multiplication, and ``g1 ** n`` means
    exponentiation of ``g1`` with the integer ``n``. ``g1 ** (-1)`` 
    is the inverse of ``g``. ``g1 / g2`` means ``g1 * g2 ** (-1)``.

    ``g1 ** g2`` means ``g2**(-1) * g1 * g2``. 

    Multiplying an object of class ``P3_node`` with an object
    of class ``AutP3`` means application of an automorphism of ``P3``
    to a point or line in ``P3``.
    """
    __slots__ = "perm"
    group_name = "Aut_P3"
    group = None       # will be set to StdAutP3Group later
    assert  find_cross([0,1,2,5])

    def __init__(self, mapping = None, data = None):
        if not mapping:
            self.perm = list(range(13))
        elif isinstance(mapping, AutP3):
            self.perm = mapping.perm[:]
        elif isinstance(mapping, str) and len(mapping) == 1:
            if mapping == 'r':
                self.perm = p3_mapping(data, random = True)
            elif  mapping in 'pl':
                self.perm = check_perm_P3([x % 13 for x in data])
                if mapping == 'l':
                    self.perm = line_map_from_map(self.perm)
            else:
                err = "Bad string in constructor of class AutP3"
                raise TypeError(err)
        else:
            self.perm = p3_mapping(mapping)

    def check(self):
        r"""Check automorphism for consistency via 'assert' statements

        ``a.check()`` returns ``a``.
        """
        check_perm_P3(self.perm)
        return self

    def __hash__(self):
        perm = self.perm
        return perm[0] + (perm[1] << 4) + (perm[2] << 8) + (perm[5] << 12)


    def order(self):
        r"""Return order of element of the group AutP3"""
        if self.perm == NEUTRAL_PERM_P3:
            return 1
        pwr = AutP3(self)
        for o in range(2, 14):
            pwr *= self
            if pwr.perm == NEUTRAL_PERM_P3:
                return o
        s = "Cannot compute order of element of AutP3"
        raise ValueError(s)

    def map(self):
        r"""Return the automorphism as a permutation list of length 26

        Element ``g`` maps P3 node``x`` to node ``g.map[x]``. Here the
        indices and values in the returned list are the numbers of the 
        nodes as in class ``P3_node``.
        """
        line_perm = [x + 13 for x in line_map_from_map(self.perm)]
        return self.perm[:] + line_perm

    def point_map(self):
        r"""Return the automorphism as a point permutation list of length 13

        Element ``g`` maps P3 point``x`` to point ``g.map[x]``. 
        """
        return self.perm[:]

    def line_map(self):
        r"""Return the automorphism as a line permutation list of length 13

        Element ``g`` maps line ``x`` to line ``g.line_map[x]``.
        The entries in the list are reduced modulo 13.
        """
        return line_map_from_map(self.perm)

    def isneutral(self):
        return self.perm == NEUTRAL_PERM_P3



@singleton
class AutP3Group(AbstractGroup):
    """Auxilary class for class ``AutP3``

    This makes the methods in class ``AbstractGroup`` available to
    instancs of class ``AutP3``.
    """
    word_type = AutP3              # type of an element (=word) in the group
    conversions = {}

    def __init__(self):
        super(AutP3Group, self).__init__()

    def __call__(*args, **kwds):
        raise TypeError("Class AutP3Group object is not callable")

    def atom(self, tag = None, data = None):
        err = "Class AutPlGroup has no attribute 'atom'"
        raise AttributeError(err)

    @staticmethod
    def _imul(g1, g2):
        return AutP3('p', mul_perm_P3(g1.perm, g2.perm))

    @staticmethod
    def _invert(g1):
        return AutP3('p', invert_perm_P3(g1.perm))

    def copy_word(self, g1):
        return AutP3(g1)

    def _equal_words(self, g1, g2):
        return g1.perm == g2.perm

    def str_word(self, g):
        """Convert group atom g to a string

        """
        return "AutP3" + str(tuple(g.perm))

StdAutP3Group = AutP3Group()   # This is the only instance of AutP3Group

AutP3.group = StdAutP3Group






         

