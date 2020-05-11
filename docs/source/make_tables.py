
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join('..', '..', 'src')))

from mmgroup import GCode, Cocode
from mmgroup.structures.gcode import str_vector24


OUTPUT_FILE = "golay_basis_data.inc"

comment = """

.. comment : This file has been created automatically 
.. comment : by executing docs/source/make_tables.py

"""


gcode_basis = """
The basis vectors ``b_0``,..., ``b_11`` of the Golay code
are given in the following table. The table also shows the
cocyles of the basis vectors as explained below::

         Basis vectors of the Golay code             Cocycle theta
         binary: bit 0,...,23            hex         bit 0,...,11
"""


cocode_basis = """

Representatives of the basis vectors ``c_0``,..., ``c_11`` of the 
Golay cocode are given in the following table::

         Basis vectors of the cocode
         binary: bit 0,...,23            hex
"""




def write_gcode_basis():
    with open(OUTPUT_FILE, "wt") as f:
        f.write(comment)
        f.write(gcode_basis)
        b = GCode.basis
        for i in range(12):
            v = str_vector24(b[i])
            n = "b_%d" % i
            h =  "0x%06x" % b[i]
            th = str_vector24(GCode(1 << i).theta().ord)[:14]
            print("  %4s:  %s,  %s;   %s" % (n, v, h, th), 
                    file = f)
        f.write(cocode_basis)
        b = Cocode.basis
        for i in range(12):
            v = str_vector24(b[i])
            n = "c_%d" % i
            h =  "0x%06x" % b[i]
            print("  %4s:  %s,  %s" % (n, v, h), file = f)
        print("", file = f)
        f.close()



if __name__ == "__main__":
    write_gcode_basis()
