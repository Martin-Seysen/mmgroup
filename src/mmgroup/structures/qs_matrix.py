

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import re
from collections.abc import Iterable
from numbers import Integral, Complex
import math
import cmath
from random import randint

import numpy as np
from mmgroup.clifford12 import QState12, as_qstate12 
from mmgroup.clifford12 import qstate12_unit_matrix 
from mmgroup.clifford12 import qstate12_column_monomial_matrix
from mmgroup.clifford12 import qstate12_row_monomial_matrix
from mmgroup.clifford12 import qstate12_pauli_matrix
from mmgroup.clifford12 import error_string
from mmgroup.clifford12 import qstate12_pauli_vector_mul
from mmgroup.clifford12 import qstate12_pauli_vector_exp

class QStateMatrix(QState12):
    """This class models a quadratic state matrix
       
    Quadratic state matrices are described in the *guide* at

    https://mmgroup.readthedocs.io/en/latest/guide.html#module-mmgroup.dev.clifford12.qstate12  .   
    
    :param rows:
    
        Binary logarithm of the number of rows of the matrix
    
    :type rows: A nonnegative ``int``
    
    :param cols:
    
        Binary logarithm of the number of columns of the matrix
    
    :type cols: A nonnegative ``int``
    
    :param data:
    
        The data of the matrix as described below
        
    :param data:
    
        Evaluated according to parameter ``data`` as described below
        
    :raise:
        * TypeError if ``type(data)`` is not as expected.
        * ValueError if ``data`` cannot be converted to an
          instance of class  ``QStateMatrix``.
    
    
    In terms of the theory of quantum computing, ``rows, cols = 0, n`` 
    creates a column vector or a *-ket* ``|v>`` corresponding to a 
    state of of ``n`` qubits, and ``rows, cols = n, 0`` creates a 
    row vector or a *-bra* ``<v|`` corresponding to a linear function
    on a state of ``n`` qubits.

    If ``rows == cols == n`` and the created ``2**n`` times ``2**n``
    matrix is invertible, then the matrix is (a scalar multiple of) 
    an element of the complex Clifford :math:`\mathcal{X}_{12}` of 
    ``n`` qubits described in :cite:`NRS01`.

    If ``rows`` is an instance of this class then a copy of 
    that instance is created.

    If ``rows`` an instance class ``QState`` then that source
    is interpreted as a *-ket* and a copy of that *-ket* is 
    created.

    If ``rows`` and ``cols`` are integers then ``data``  may be:

      * ``None`` (default). Then the zero matrix is created.

      * An integer ``v``, if ``rows`` or ``cols`` is ``0``. 
        Then the state is set to ``|v>`` or ``<v|``,  respectively.

      * The integer ``1`` if ``rows`` == ``cols``. Then a unit 
        matrix is created.

      * A list of integers. Then that list of integers must encode 
        a valid pair ``(A, Q)`` of bit matrices that make up a 
        state as in class ``QState``, as dscribed in the *guide*. 
        In this case parameter ``mode`` is evaluated as follows:
              
          * 1: create matrix ``Q`` from lower triangular part
              
          * 2: create matrix ``Q`` from upper triangular part
               
          * Anything else: matrix ``Q`` must be symmetric.
          
    As in ``numpy``, matrix multiplication of quadratic state 
    matrices is done with the ``@`` operator and elementwise
    multiplication of such matrices is done with the ``*``
    operator. A quadratic state matrix may also be multiplied 
    by a scalar. Here the scalar must be zero or of the form

    .. math::
         2^{e/2} \cdot w \mid e \in \mathbb{Z}, \;
         w \in \mathbb{C}, \, w^8 = 1   \; .
   
    A matrix of type ``QStateMatrix`` may be indexed with square
    brackets as in ``numpy`` in order to obtain entries, rows, 
    columns or submatrices. Then a complex ``numpy`` array (or a 
    complex number) is returned as in ``numpy``. It is  not 
    possible to change the matrix via item assignment. So the 
    easiest way to obtain a complex version of an instance ``m`` 
    of type `QStateMatrix`` is to write ``m[:,:]``.
    
    Officially, we support matrices with ``rows, cols <= 12``
    only. Methods of this class might work for slightly 
    larger matrices. Any attempt to constuct a too large
    matrix raises ValueError.
    """
    UNDEF_ROW = 255
    
    def __init__(self, rows, cols = None, data = None, mode = 0):
        if isinstance(cols, Integral):
            if isinstance(rows, Integral):
                n = rows + cols
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
                self._reshape(rows, cols)
        elif cols is None:
            source = rows
            if isinstance(source, QState12):
                super(QStateMatrix, self).__init__(as_qstate12(source))
            elif isinstance(source, QStateMatrix):
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
    


    def reshape(self, shape = (), copy = True):
        """Reshape matrix to given ``shape``

        ``shape[0] + shape[1] = self.shape[0] + self.shape[0]`` 
        must hold.

        If one of the values ``shape[0]``, ``shape[1]`` is negative,
        the other value is calculated from the sum 
        ``self.shape[0] + self.shape[1]``. 

        Shape defaults to ``(0, -1)``.
        """  
        m = QStateMatrix(self) if copy else self
        if isinstance(shape, int):
            shape = (shape, -1)
        while len(shape) < 2:
            shape += (-1),
        m._reshape(shape[0], shape[1])
        return m
        
    def conjugate(self, copy = True):
        """Return complex conjugate of the matrix
        
        Returns an instance of class ``QStateMatrix``.
        """
        m = QStateMatrix(self) if copy else self   
        return super(QStateMatrix, m).conjugate()   

    conj = conjugate        

       
    def rot_bits(self, rot, nrot, n0 = 0, copy = True):    
        """Wrapper for corresponding function in class QState12"""
        m = QStateMatrix(self) if copy else self
        super(QStateMatrix, m).rot_bits(rot, nrot, n0)   
        return m   
        
    def xch_bits(self, sh, mask, copy = True):    
        """Wrapper for corresponding function in class QState12"""
        m = QStateMatrix(self) if copy else self
        super(QStateMatrix, m).xch_bits(sh, mask)   
        return m   
        
    def extend(self, j,nqb, copy = True):    
        """Wrapper for corresponding function in class QState12
        
        The function returns *ket*, i.e. a column vector. 
        """
        m = QStateMatrix(self) if copy else self
        super(QStateMatrix, m).extend(j, nqb) 
        m._reshape(0, m.ncols)       
        return m   

    def extend_zero(self, j, nqb, copy = True):    
        """Wrapper for corresponding function in class QState12
        
        The function returns *ket*, i.e. a column vector. 
        """
        m = QStateMatrix(self) if copy else self
        super(QStateMatrix, m).extend_zero(j, nqb) 
        m._reshape(0, m.ncols)        
        return m   

    def restrict(self, j,nqb, copy = True):    
        """Wrapper for corresponding function in class QState12
        
        The function returns *ket*, i.e. a column vector. 
        """
        m = QStateMatrix(self) if copy else self
        super(QStateMatrix, m).restrict(j, nqb) 
        m._reshape(0, m.ncols)     
        return m   

    def restrict_zero(self, j, nqb, copy = True):    
        """Wrapper for corresponding function in class QState12
        
        The function returns *ket*, i.e. a column vector. 
        """
        m = QStateMatrix(self) if copy else self
        super(QStateMatrix, m).restrict_zero(j, nqb) 
        m._reshape(0, m.ncols)       
        return m   

    def sumup(self, j,nqb, copy = True):    
        """Wrapper for corresponding function in class QState12
        
        The function returns *ket*, i.e. a column vector. 
        """
        m = QStateMatrix(self) if copy else self
        super(QStateMatrix, m).sumup(j, nqb) 
        m._reshape(0, m.ncols)        
        return m   


    @property 
    def T(self):
        """Return transposed matrix as in numpy"""   
        return self.copy().transpose()

    @property 
    def H(self):
        """Return conjugate transposed matrix as in numpy"""   
        m = self.copy().transpose()
        return super(QStateMatrix, m).conjugate()    
     
    def complex(self):
        """Return complex matrix as a numpy array
        
        ``m.complex()`` is equivalent to ```m[:,:]``.
        """
        a = super(QStateMatrix, self).complex()
        a = a.reshape((1 << self.shape[0], 1 << self.shape[1])) 
        return a        

    def complex_unreduced(self):
        """For tests only: Return complex matrix of state 
        
        Returns same result as method reduce().

        The standard method reduce() creates a reduced copy before 
        calculating the complex matrix, which is usually much 
        faster. This method does not reduce the matrix.
        """
        a = super(QStateMatrix, self).complex_unreduced()
        a = a.reshape((1 << self.shape[0], 1 << self.shape[1])) 
        return a 

       
    def lb_norm2(self):    
        """Return binary logarithm of squared operator norm.
        
        The operator norm is the largest absolute singular value
        of the matrix.
        
        Return -1 if matrix is zero.
        """
        m = self.copy().reduce()
        if m.nrows == 0:
            return -1
        return m.factor[0] + m.nrows - 1 - m.lb_rank() 

    def inv(self):
        """Return inverse matrix
     
        Returns an instance of class ``QStateMatrix``.
        Raise ValueError if matrix is not invertible.
        """
        return self.copy()._mat_inv()
        

    def pauli_vector(self):
        try:
            return super(QStateMatrix,self).pauli_vector()
        except ValueError:
            err = "QStateMatrix object is not in the Pauli group"
            print("\n%s:\n%s\n" % (err, str(self)))
            raise            
            
    def __matmul__(self, other):
        p1, p2 = self.copy(), other.copy()
        try:
            return p1._matmul(p2)
        except ValueError:
            err = "Multiplying QStateMatrix objects of shape %s and %s"
            if isinstance(other, QStateMatrix):
                print("\n" + err % (self.shape, other.shape))
            raise          

        
    
    def  __mul__(self, other):
        if isinstance(other, Complex):
            e, phi = complex_to_qs_scalar(other)
            if e is None:
                return self.copy().set_zero()
            else:
                return self.copy().mul_scalar(e, phi)
        elif isinstance(other, QStateMatrix):
            if self.shape == other.shape:
                p1, p2 = self.copy(), other.copy()
                p1._qstate12_product(p2, self.ncols, 0)
                p1._reshape(*self.shape)
                return p1
            else:
                err = "QStateMatrix instances must have same shape"
                raise ValueError
        else:
            err = "Byte type for multiplication with QStateMatrix"
            raise TypeError(err)
            
    __rmul__ = __mul__       

    def __truediv__(self, other):
        if isinstance(other, Complex):
            return self.__mul__(1.0/other)
        else:
            err = "Byte type for multiplication with QStateMatrix"
            raise TypeError(err)
            
    def __neg__(self):
        return self.copy().mul_scalar(0, 4)    

    def __pos__(self):
        return self    
        
    def __eq__(self, other):
        try:
            return self._equal(other)
        except ValueError:
            if isinstance(other, QStateMatrix):
                err = "Comparing QStateMatrix objects of shape %s and %s"
                print("\n" + err % (self.shape, other.shape))
            raise          
        
    def __getitem__(self, item):
        if not isinstance(item, tuple):
            item = (item,)
        while len(item) < 2:
            item = item + (None,)
        a0, f0 = _as_index_array(item[0], self.shape[0]) 
        a1, f1 = _as_index_array(item[1], self.shape[1])
        if f0 & f1:
            return self.complex()
        a0 = a0 << self.shape[1]
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

def rand_qs_matrix(rows, cols, data_rows):
    limit = (1 << (rows  + cols + data_rows)) - 1 
    data = [randint(0, limit) for i in range(data_rows)]
    return QStateMatrix(rows, cols, data, mode = 1)

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
        
    
def format_state(q):
    """Return  a ``QStateMatrix`` object as a string."""
    rows, cols = q.shape
    try:
        data = q.data
    except ValueError:
        print("\nError: Bad instance of class QStateMatrix:\n") 
        print(format_state_raw(q))        
        raise    
        
    e = q.factor
    str_e = format_scalar(*e) if len (data) else "0"
    str_data = format_data(data, rows, cols, reduced = False)   
    qtype = STATE_TYPE[bool(rows), bool(cols)]
    s = "<%s %s" % (qtype, str_e)
    if len(str_data):
       s += " *\n" + str_data 
    s += ">\n"                  
    return s

def format_state_raw(q):
    s = ("QStateMatrix, shape = %s, factor = %s, nrows = %s, cols = %s\n"
        % (q.shape, q.factor, q.nrows, q.ncols))
    s += "data:\n"
    length = q.nrows + q.ncols
    for d in q.raw_data[:q.nrows]:
        s += binary(d, 0, length) + "\n"
    return s
        

####################################################################
# Computing an array of indices
####################################################################
        
        
def _as_index_array(data, nqb):
    """Convert an index ``data`` to an array of indices.
    
    ``data`` in an object for indexing a one-dimesnional numpy
    array of length ``nqb``. 
    
    Return a pair ``(a, full)`` where ``a`` is a numpy array 
    containing the indices and ``full`` is True iff
    ``a`` contains the data ``range(1 << nqb)``.
    """
    mask = (1 << nqb) - 1
    if isinstance(data, Integral):
        return np.array(data & mask, dtype = np.uint32), False
    if data is None:
        return np.arange(1 << nqb, dtype = np.uint32), True
    if isinstance(data, slice):
        ind = data.indices(1 << nqb)
        full = ind == (0, 1 << nqb, 1) 
        return np.arange(*ind, dtype = np.uint32), full
    ind = np.array(data, dype = np.uint32, copy = False) 
    if len(ind.shape) > 1:
        err = "Bad index type for QState12 array"
        raise TypeError(err)
    return  ind & mask, False 
    
    
    
####################################################################
# convert a complex number to a factor for class QStateMatrix 
####################################################################
    

MIN_ABS = 2.0**-1024
EPS = 1.0e-8 
    
def complex_to_qs_scalar(x):
    """Convert complex number to scalar factor
    
    Return ``(e, phi)`` if ``x == 2**(0.5*e) * z**phi``
    with ``z = (1+1j)/sqrt(2)``.
    
    Return ``None, None`` if ``x == 0``
    
    raise ValueError otherwise. We accept a relative error of
    about ``1.0e-8`` for ``x``.
    
    So it is safe to write e.g. ``2**-0.5 * (-1+1j) * m`` 
    if ``m`` is an  instance of class ``QStateMatrix``.    
    """
    r, phi = cmath.polar(x)
    if r <= MIN_ABS:
        return None, None
    e = 2 * (math.log(r) / math.log(2))
    phi8 = 4 * phi / math.pi
    e_r, phi8_r = round(e), round(phi8)
    if (max(abs(e - e_r), abs(phi8 - phi8_r)) < EPS):
        return e_r, phi8_r
    err = "Cannot convert number to factor for QStateMatrix"
    raise ValueError(err)    

####################################################################
# Some wrappers
####################################################################

    
def flat_product(a, b, nqb, nc):
    return a.copy()._qstate12_product(b.copy(), nqb, nc)
    
        
def qstate_column_monomial_matrix(data):
    nqb = len(data) - 1
    qs = QStateMatrix(nqb, nqb, 1) 
    qstate12_column_monomial_matrix(qs, nqb, data)
    return qs
    
def qstate_row_monomial_matrix(data):
    nqb = len(data) - 1
    qs = QStateMatrix(nqb, nqb) 
    qstate12_row_monomial_matrix(qs, nqb, data)
    return qs

def qstate_unit_matrix(nqb):
    qs = QStateMatrix(nqb, nqb) 
    qstate12_unit_matrix(qs, nqb)
    return qs

def qstate_pauli_matrix(nqb, v):
    qs = QStateMatrix(nqb, nqb)
    qstate12_pauli_matrix(qs, nqb, v)    
    return qs
    
def qstate_ctrl_not_matrix(nqb, vc, v):
    """Return Transformation matrix for ctrl-not gate
    
    ``nqb`` is the rank of the matrix.
    TODO: ``vc, v`` yet to be documented.

    """
    qs = qstate_unit_matrix(nqb)
    mask = (1 << nqb) - 1
    qs.gate_ctrl_not(vc & mask, v & mask) 
    return qs.reduce()    

def qstate_phi_matrix(nqb, v, phi):
    qs = qstate_unit_matrix(nqb)
    qs.gate_phi(v << nqb, phi) 
    return qs    

def qstate_ctrl_phi_matrix(nqb, v1, v2):
    qs = qstate_unit_matrix(nqb)
    qs.gate_ctrl_phi(v1 << nqb, v2) 
    return qs    

def qstate_hadamard_matrix(nqb, v):
    qs = QStateMatrix(nqb, nqb, 1)
    qs.gate_h(v)    
    return qs.reduce()

def qstate_pauli_vector_mul(nqb, v1, v2):
    return qstate12_pauli_vector_mul(nqb, v1, v2)

def qstate_pauli_vector_exp(nqb, v, e):
    return qstate12_pauli_vector_exp(nqb, v, e)
