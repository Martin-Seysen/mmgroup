

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import re
from collections.abc import Iterable
from numbers import Integral

import numpy as np
from mmgroup.clifford12 import QState12 




class QStateMatrix(QState12):
    def __init__(self, rows, cols = None, data = None):
        """Create a 2**rows times 2**cols quadratic state matrix

        If ``rows`` and ``cols`` are integers then ``data``
        may be:

            * ``None`` (default). Then the zero matrix is created.

            * An integer ``v``, if ``rows`` or ``cols`` is ``0``. 
              Then the state is set to ``|v>`` or ``<v|``, 
              respectively.

            * The integer ``1`` if ``rows`` == ``cols``. Then 
              a unit matrix is created.

            * A list of integers. Then that list of integers must 
              encode a valid pair ``(A, Q)`` of bit matrices that 
              make up a state as in class ``QState``. 

        So ``rows, cols = 0, n`` creates a column vector or a *-ket*
        ``|v>`` corresponding to a state of of ``n`` qubits, and 
        ``rows, cols = n, 0`` creates a row vector or a *-bra* ``<v|`` 
        corresponding to a linear function on a state of ``n`` qubits.

        If ``source`` is an instance of this class then a copy of 
        that instance is created.

        If ``source`` an instance class ``QState`` then that source
        is interpreted as a *-ket* and a copy of that *-ket* is 
        created.
        """
        if isinstance(cols, Integral):
            if isinstance(rows, Integral):
                self.rows, self.cols, n = rows, cols, rows + cols
                if data is None or isinstance(data, Iterable):
                    super(QStateMatrix, self).__init__(n, data)
                elif isinstance(data, Integral):
                    if rows * cols == 0:
                        super(QStateMatrix, self).__init__(n, data)
                    elif data == 1:
                        err = "Unit matrix not implemented"
                        raise ValueError(err)
                else:
                    err = "Bad data type for  QStateMatrix"
                    raise TypeError(err) 
        elif cols is None:
            source = rows
            if isinstance(source, QStateMatrix):
                self.rows, self.cols = source.rows, source.cols
                super(QStateMatrix, self).__init__(source)
            elif isinstance(source, QState):
                self.rows, self.cols = 0, source.ncols
                super(QStateMatrix, self).__init__(source)
            else:
                err = "Illegal source type for QStateMatrix"
                raise TypeError(err) 
        else:
            err = "Cannot construct QStateMatrix from given objects" 
            raise TypeError(err) 
                
              

    def __str__(self):
        return format_state(self)


####################################################################
# Formatting a QStateMatrix object
####################################################################


PHASE = [ ("",""),  (""," * (1+1j)"),  (""," * 1j"), ("-"," * (1-1j)"),
          ("-",""), ("-"," * (1+1j)"), ("-"," * 1j"), (""," * (1-1j)"), 
]

def format_scalar(e, phase):
    if (phase & 1): e -= 1
    prefix, suffix = PHASE[phase]
    if 0 <= e <= 12 and not e & 1:
        s_exp = str(int(2**(e/2)))
    else:
        s_exp = ("2**%.1f" if e & 1 else "2**%d") % (e/2)
        if prefix:
            s_exp = "(" + s_exp + ")"
    return prefix + s_exp + suffix  
  
BRACKETS = { 
   (0,0):("  <","scalar",">"),   # a scalar
   (0,1):("  <","","|"),         # a row vector or a  *bra-*
   (1,0):("  |","",">"),         # a column vector or a  *-ket*
   (1,1):("  |", "><","|")       # a matrix
}


def binary(n, start, length, leading_zeros = True):
    if length == 0:
        return ""
    n = (n >> start) & ((1 << length) - 1)
    b = format(n, "b")
    c = "0" if leading_zeros else " "
    return c * (length - len(b)) + b

def binary_q(n, start, length, pos):
    def pm_repl(m):
        return "-" if  m.group(0) == "1" else "+"
    def j_repl(m):
        return "j" if  m.group(0) == "1" else "."
    if length < 2:
        return ""
    s =  binary(n, start, length) 
    s = "".join(s[::-1])
    if (pos):
       s = re.sub("[01]", pm_repl, s, pos)
    if (pos < len(s)):
       s = re.sub("[01]", j_repl, s, 1) 
    s = re.sub("[01]", pm_repl, s)
    return s

                      
def format_data(data, rows, cols, reduced = False):
    s = ""
    nrows = len(data)
    if len(data) < 2 and  rows + cols == 0:
        return s
    if len(data) == 0:
         data = [0]
    left, mid, right = BRACKETS[bool(rows), bool(cols)]
    left_bl = " " * len(left)
    mid_bl = " " * len(mid)
    right_bl = " " * len(right)
    for i, d in enumerate(data):
        c = binary(d, 0, cols, not reduced) 
        r = binary(d, cols, rows, not reduced) 
        q = binary_q(d, rows + cols, len(data), i)
        s += left + r + mid + c + right + " " + q + "\n"
        left, mid, right = left_bl, mid_bl, right_bl
    return s
        

STATE_TYPE = { (0,0) : ("QState scalar"), 
               (0,1) : ("QState row vector"), 
               (1,0) : ("QState column vector"),
               (1,1) : ("QState matrix")
}    
        
    
def format_state(q, reduced = False):
    q.check()
    rows, cols = q.rows, q.cols
    data = q.data
    e = q.factor
    str_e = format_scalar(*e) if len (data) else "0"
    str_data = format_data(data, rows, cols, reduced = False)   
    qtype = STATE_TYPE[bool(rows), bool(cols)]
    s = "<%s %s" % (qtype, str_e)
    if len(str_data):
       s += " *\n" + str_data 
    return s + ">\n"                  
        