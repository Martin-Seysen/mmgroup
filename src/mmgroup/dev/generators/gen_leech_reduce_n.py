r"""This is yet a stub


"""
from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

from mmgroup.bitfunctions import bw24, lmap, pivot_binary_high
from operator import __or__, __xor__

import random
import numpy
from functools import reduce
from numpy import array, zeros, uint8, uint16, uint32

from mmgroup.dev.mat24.mat24_ref import Mat24


#######################################################################
# Auxiliary functions 
#######################################################################



def list_to_vect(list_):
    return  reduce(__or__, [1 << x for x in list_], 0) 


def y_table(g, cocodes, sign, code_omega, cocode_x, verbose = 0):
    """Return table for function  ``gen_leech2_reduce_n`` 


    This table is used to map an entry :math:`x_d x_\delta`,
    where :math:`d` is a (signed) octad or dodecad, and 
    where :math:`\delta` is even, to a standard form. The
    table is to be used for a fixed Golay code element 
    :math:`d` up to sign and up to a factor :math:`\Omega`.

    Such a mapping is done via conjugation with an element
    of shape  :math:`y_e x_\epsilon`, with math:`e` a Goly
    code word math:`\epsilon` an even cocode word..

    Here :math:`d` is given by parameter ``g``. Thee we choose
    a Golay code word ``e`` with :math:`g \cap e = \delta`. 
    Using the returned table, ``e`` may be calculated as

    :math:`e = \sum_{i=0}^10 \delta_i \cdot t_i` ,

    where :math:`\delta_i` is the ``i``-th bit of
    :math:`\delta`, and :math:`t_i` is the entry ``table[i+1]`` 
    of the table.
    
    ``table[0]`` is equal to :math:`\theta(d)`. ``table[12]`` 
    is a sign bit ``s``. If the bit 23 of :math:`v \cdot y_e`
    differs from ``s`` then we have to replace :math:`e` by
    :math:`e + e'`, where :math:`e'` is given by ``table[13]``. 
    This operation adds  :math:`x_\Omega` to  :math:`x_d` in
    the cases relevant for function ``gen_leech2_reduce_n``.
    If :math:`v \cdot y_{e + e'}` has a sign bit set then
    we have to multiply that vector with  :math:`x_\epsilon`,
    where the cocode word :math:`epsilon` is given by 
    ``table[13]``. Entries ``table[11,12,13]`` are copied
    from parameters ``sign, code_omega, cocode_x``,
    respectively.

    Parameter ``cocode_x`` is an list of cocode words
    that will be meapped to zero by applying the table.

    On input, all Golay code and cocode word must be given
    as list of integers representing the bit positions set.

    In the returned table, all Golay code or cocode words
    are given as bit vector in **Golay code** or **cocode**
    representation, respectively.
    """
    gv = Mat24.vect_to_gcode(list_to_vect(g))
    b = []
    for i in range(11):
        b.append((Mat24.ploop_cap(gv, 1 << i) << 12) + (1 << i))
    for c in cocodes:
        b.append(Mat24.vect_to_cocode(list_to_vect(c)) << 12)
    if verbose:
        print("Function y_table, g =", hex(list_to_vect(g)))
        print("cocodes =", cocodes)
        print("cocode_x =", cocode_x)
        print("b =", [hex(x) for x in b])
    b, columns = pivot_binary_high(b)
    if verbose:
        print("reduced b =", [hex(x) for x in b])
        print("columns =", columns)
    table = [Mat24.ploop_theta(gv)]
    table1 = []
    for i in range(10, -1, -1):
        if columns[0] == i + 12:
            table1.append(b[0] & 0xfff)
            b, columns = b[1:], columns[1:]
        else:
            table1.append(0)
    table1.reverse()
    table += table1
    table.append(sign)
    table.append(Mat24.vect_to_gcode(list_to_vect(code_omega)))
    table.append(Mat24.vect_to_cocode(list_to_vect(cocode_x)))
    if verbose:
        print("table =", [hex(x) for x in table])
    return table



#######################################################################
# GenLeechReduceY
#######################################################################



class GenLeechReduceY(object):
    """This is yet a stub

    """

    tables = {}

    directives = {}


