import sys
import os


if __name__ == "__main__":
    sys.path.append("../../../")

from mmgroup.dev.mat24.mat24_ref import Mat24 as mat24
from mmgroup.dev.mat24.make_addition_table import make_addition_operations
from mmgroup.dev.mat24.make_addition_table import check_addition_operations
from mmgroup.dev.mat24.make_addition_table import display_addition_operations

def make_short_ops_table(verbose = 0):
    if verbose:
        s = "Map cocode basis to short cocode vectors [0,i], 1 <= i < 24"
        print(s)
    m = [0] + [mat24.vect_to_cocode(1 + (1 << i)) for i in range(1,24)]
    for i in range(1,24):
        cc = mat24.cocode_syndrome(m[i], 0)
        assert cc == (1 << i) + 1
        if verbose:
            print("0x%03x : 0x%06x" % (m[i], cc))
    singleton = True
    ops = make_addition_operations(m, singleton, granularity = 4)
    check_addition_operations(m, ops, singleton)
    if verbose:
        display_addition_operations(ops)
        print("%d operations" % len([x for x in ops if len(x) >= 4]))
    return ops


def C_table_from_short_ops(ops):
    in_table = []  
    op_table = []
    regs = {}
    for op in ops:
        if len(op) < 4:
            if op[1] > 0:
                assert op[1] == 1 << len(in_table)
                in_table.append(op[0])
                regs[op[1]] = op[0]
        else:
            op_table += [regs[x] for x in op[2:]]
            op_table.append(op[0])
            regs[op[1]] = op[0]
    assert max(op_table) < 24
    assert len(op_table) % 3 == 0
    assert len(ops) == len(in_table) + len(op_table) // 3 + 1
    return in_table, op_table


class ShortCocodeTables:
    ops = make_short_ops_table(verbose = 0)
    in_table, op_table = C_table_from_short_ops(ops)

    tables = {
       "ShortCocode_InTable" : in_table,
       "ShortCocode_OpTable" : op_table,
       "ShortCocode_LenOpTable" : len(op_table),
       "ShortCocode_NumOps" : len(op_table) // 3,
    }

    directives = {}


    @staticmethod
    def display_op():
        make_short_ops_table(verbose = 1)


class Tables(ShortCocodeTables):
    pass



if __name__ == "__main__":
    ShortCocodeTables.display_op()
    print("Input table:")
    print(ShortCocodeTables.in_table)
    print("Operation table:")
    print(ShortCocodeTables.op_table)
   




