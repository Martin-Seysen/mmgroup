

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import re
from collections.abc import Iterable
from numbers import Integral
from random import randint

import numpy as np
from mmgroup.clifford12 import QState12, as_qstate12 
from mmgroup.clifford12 import qstate12_unit_matrix 
from mmgroup.clifford12 import qstate12_matmul, qstate12_prep_mul
from mmgroup.clifford12 import qstate12_product


class QStateMatrix(QState12):
    def __init__(self, rows, cols = None, data = None, mode = 0):
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
              make up a state as in class ``QState``. In this case 
              parameter ``mode`` is evaluated as follows:
              
               * 1: create matrix ``Q`` from lower triangular part
              
               * 2: create matrix ``Q`` from upper triangular part
               
               * Anything else: matrix ``Q`` must be symmetric.

        So ``rows, cols = 0, n`` creates a column vector or a *-ket*
        ``|v>`` corresponding to a state of of ``n`` qubits, and 
        ``rows, cols = n, 0`` creates a row vector or a *-bra* ``<v|`` 
        corresponding to a linear function on a state of ``n`` qubits.

        If ``rows`` is an instance of this class then a copy of 
        that instance is created.

        If ``rows`` an instance class ``QState`` then that source
        is interpreted as a *-ket* and a copy of that *-ket* is 
        created.
        """
        if isinstance(cols, Integral):
            if isinstance(rows, Integral):
                self.rows, self.cols, n = rows, cols, rows + cols
                if data is None or isinstance(data, Iterable):
                    super(QStateMatrix, self).__init__(n, data, mode)
                elif isinstance(data, Integral):
                    if rows * cols == 0:
                        super(QStateMatrix, self).__init__(n, data)
                    elif rows == cols and data == 1:
                        qstate12_unit_matrix(self, rows)
                else:
                    err = "Bad data type for  QStateMatrix"
                    raise TypeError(err) 
        elif cols is None:
            source = rows
            if isinstance(source, QStateMatrix):
                self.rows, self.cols = source.rows, source.cols
                super(QStateMatrix, self).__init__(as_qstate12(source))
            elif isinstance(source, QState12):
                self.rows, self.cols = 0, source.ncols
                super(QStateMatrix, self).__init__(as_qstate12(source))
            else:
                err = "Illegal source type for QStateMatrix"
                raise TypeError(err) 
        else:
            err = "Cannot construct QStateMatrix from given objects" 
            raise TypeError(err) 
            
    def copy(self):
        """Return a copy of the matrix"""    
        return QStateMatrix(self)   
    
    @property
    def shape(self):
        return (self.rows, self.cols)  


    def reshape(self, shape = (), copy = True):
        """Reshape matrix to given ``shape``

        ``shape[0] + shape[1] = self.rows + self.cols`` must hold.

        If on of the values ``shape[0]``, ``shape[1]`` is negative,
        the other value is calculated from ``self.rows + self.cols``. 

        Shape default to ``(-1, 0)``.
        """  
        m = QStateMatrix(self) if copy else self
        if isinstance(shape, Integral): 
            shape = (shape,)
        while len(shape) < 2:
            shape += (-1 if min(shape) >= 0 else 0),
        rows, cols = shape
        if rows < 0:
            rows = m.ncols - cols
        if cols < 0:
            cols = m.ncols - rows
        if rows + cols != m.ncols or min(rows, cols) < 0:
            err = "Bad shape for reshaping  QStateMatrix"
            raise ValueError(err)
        m.rows, m.cols = rows, cols
        return m
     
    def complex(self):
        """Return complex matrix of state as numpy array"""
        a = super(QStateMatrix, self).complex()
        a = a.reshape((1 << self.rows, 1 << self.cols)) 
        return a        

    def complex_unreduced(self):
        """For tests only: Return complex matrix of state 
        
        Returns same result as method reduce().

        The standard method reduce() creates a reduced copy before 
        calculating the complex matrix, which is usually much 
        faster. This method does not reduce the matrix.
        """
        a = super(QStateMatrix, self).complex_unreduced()
        a = a.reshape((1 << self.rows, 1 << self.cols)) 
        return a  

    def __matmul__(self, other):
        r1, c1 = self.shape
        r2, c2 = other.shape
        if r2 != c1:
            err = "Shape mismatch in QStateMatrix multiplication"
            raise ValueError(err)
        result = self.copy()
        qstate12_matmul(result, other.copy(), c1)
        result = QStateMatrix(result)
        result.reshape((r1, c2), copy = False)
        return result
            
        
    def __getitem__(self, item):
        if not isinstance(item, tuple):
            item = (item,)
        while len(item) < 2:
            item = item + (None,)
        a0 = _as_index_array(item[0], self.shape[0]) << self.shape[1]
        a1 = _as_index_array(item[1], self.shape[1])
        shape_ =  a0.shape + a1.shape 
        a = np.ravel(a0)[:, np.newaxis] + np.ravel(a1)[ np.newaxis, :] 
        if a.dtype != np.uint32:
            a = np.array(a, dtype = np.uint32)
        a = np.ravel(a, order = "C")
        c = self.entries(a)  
        if len(shape_):        
            return c.reshape(shape_)
        return c[0]

    def __str__(self):
        return format_state(self)



####################################################################
# Creating a random a QStateMatrix object
####################################################################

def rand_qs_matrix(row, cols, data_rows):
    limit = (1 << (row  + cols + data_rows)) - 1 
    data = [randint(0, limit) for i in range(data_rows)]
    return QStateMatrix(row, cols, data, mode = 1)

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
        
        

####################################################################
# Computing an array of indices
####################################################################
        
        
def _as_index_array(data, nqb):
    mask = (1 << nqb) - 1
    if isinstance(data, Integral):
        return np.array(data & mask, dtype = np.uint32)
    if data is None:
        return np.arange(1 << nqb, dtype = np.uint32)
    if isinstance(data, slice):
        return np.arange(*data.indices(1 << nqb), 
            dtype = np.uint32)
    ind = np.array(data, dype = np.uint32, copy = False) 
    if len(ind.shape) > 1:
        err = "Bad index type for QState12 array"
        raise TypeError(err)
    return  ind & mask  

####################################################################
# Some wrappers
####################################################################

def prep_mul(a, b, nqb = None):
    if nqb is None and a.cols == b.cols:
        nqb = a.cols
    a, b = QState12(a), QState12(b)
    qstate12_prep_mul(a, b, nqb)
    return QStateMatrix(a), QStateMatrix(b)
    
    
def flat_product(a, b, nqb, nc):
    a, b = QState12(a), QState12(b)
    qstate12_product(a, b, nqb, nc)
    return QStateMatrix(a)  



    