r"""Generate table for Griess algebra multuplication by standard axis

This module generates a table for multiplying a vector in the
representation \rho by the standard axis v^+. Here v^+ is the axis
corresponding to the Golay code covector \beta = [2,3], considered
as an involution in the subgroup Q_x0 of the Monster. The operation
of v^+ on \rho is well known. The motst important case where a table
is helpful is the subspace of \rho labelled by 'T'. The entries are
labelled by a pair (o, d), where o is an octad and d an even subset
of o. Here the following cases may happen:

Case 1: Scalar product <o, \beta> = 1

   Then \beta * <o,d> = <o,d>   (assuming suitable scaling)

Case 2: Scalar product <o, \beta> = 0, d not a subset of o

    Then \beta * <o,d> = 0

Case 3: Scalar product <o, \beta> = 0, d is a subset of o

    Then <\beta * (o, d - \beta)> = <8 * (o, d - \beta)>
    if |d \cap \beta| is odd) and this is zero otherwise.
    Also, \beta * (o, d + \beta) = 0 in all cases. 

Clearly, case 3 is the most complicated case. Here it comes handy
that the mapping d \mapsto |d \cap \beta| is linear in d.

Function mat24.octad_entries(o) returns the 8 (not neccessarly sorted)
list of entries of octad o used for the encoding of part 'T' of a
vector in \rho. Here it is helpful that the two entries 2 an 3
(if present) occur in that list either at positions 0 and 1, or
at position 2 and 3, in any order.

So we may precompute the positions of the suboctads d mapping to
nonzero entries in these two cases.
"""

if __name__ == "__main__":
    import os, sys
    script_path = os.path.dirname(os.path.abspath(__file__))
    sys.path.append(os.path.join(script_path, '..', '..', '..'))
    from mmgroup import mat24

import numpy as np


###########################################################################
# 
###########################################################################

NUM_OCTADS = 759
STD_COCODE = 0x200

"""
By conincidence, if both, entries 2 and 3 are in a suboctad then
these two entries are either in positions 0 and 1 or in positions
2 and 3, in an unspecified order.
"""

INDEX_OF_3_GIVEN_INDEX_OF_2 = [1,0,3,2,None,None,None,None]

def octad_to_std_axis_op(o):
    v = mat24.octad_to_gcode(o)
    assert 0 <= v < 0xfff
    if v & STD_COCODE:
         return 1
    entries = mat24.octad_entries(o)
    if 2 in entries:
        ind2 = entries.index(2)
        pos3 = INDEX_OF_3_GIVEN_INDEX_OF_2[ind2]
        assert pos3 is not None and  entries[pos3] == 3, (ind2, pos3, entries)
        return 2 + (ind2 >= 2)
    else:
        assert 3 not in entries
        return 0


def make_table_octad_to_std_axis_op():
    assert STD_COCODE == mat24.vect_to_cocode((1 << 2) + (1 << 3))
    a = np.zeros((NUM_OCTADS + 31) // 32, dtype = np.uint64)
    for o in range(NUM_OCTADS):
        c = octad_to_std_axis_op(o)
        #assert 0 <= c < 4
        a_new = int (a[o >> 5]) + (c << (2 * (o & 0x1f)))
        a[o >> 5] = a_new  # numnpy is a bit stupid here (?)
    return a

######################################################################
# Summarizing the tables given above
######################################################################

class Tables:
    directives = {}
    def __init__(self, **kwds):
        global mat24
        from mmgroup import mat24
        self.p = p = int(kwds.get('p', 3))
        self.tables = {
            "TABLE_OCTAD_TO_STD_AX_OP" : make_table_octad_to_std_axis_op()
        }

class MockupTables:
    tables = {
        "TABLE_OCTAD_TO_STD_AX_OP" : [0]
    }
    directives = {}
    def __init__(self, **kwds):
        pass

######################################################################
# Test functions
######################################################################



if __name__ == "__main__":
    for x in Tables().tables["TABLE_OCTAD_TO_STD_AX_OP"]:
        print("0x%016x" % int(x))
    T = Tables().tables["TABLE_OCTAD_TO_STD_AX_OP"]
    for o in range(NUM_OCTADS):
        value = (int(T[o >> 5]) >> (2 * (o & 31))) & 3
        ref =  octad_to_std_axis_op(o) 
        assert value == ref, (o, value, ref) 
    1/0



