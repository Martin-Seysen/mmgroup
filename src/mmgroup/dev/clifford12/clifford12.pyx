# cython: language_level=3
"""Simulation of the  Clifford group of up to 12 qubits.

"""

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

from collections.abc import Iterable
from numbers import Integral, Complex
import cmath

import numpy as np
from libc.stdint cimport uint64_t, uint32_t, int32_t, uint8_t, int8_t
from libc.math cimport fabs, log, round


cimport cython
include "../pxd_files/clifford12.pxi"
cimport clifford12 as cl
ctypedef qstate12_type *p_qstate12_type
 
 
####################################################################
# Handling errors from C functions in file qstate.c
####################################################################
 

# Messages for error codes returned by the functions in module qstate12.c
QSTATE12_ERROR_STRINGS = {
 -1: "No bit with requested property found in QStateMatrix",
 -2: "Inconststent data in QStateMatrix",
 -3: "Qubit index error at QStateMatrix",
 -4: "State or matrix is too large for a QStateMatrix",
 -5: "Internal buffer overflow in a QStateMatrix",
 -6: "Bit matrix part Q of QStateMatrix is not symmetric",
 -7: "Bad row index for bit matrix in a QStateMatrix",
 -8: "Internal parameter error in a QStateMatrix",
 -9: "Overflow in scalar factor in a QStateMatrix",
 -10: "A qubit in a ctrl-not gate cannot control itself",
 -11: "Shape mismatch in QStateMatrix comparison",
 -101: "Shape mismatch in QStateMatrix operation",
 -102: "QStateMatrix is not invertible",
 -103: "QStateMatrix is not in the Pauli group",
 -104: "QStateMatrix is not monomial",
 
 # Additional error codes taken from file xsp2co1.c
 -201: "Internal error in operation of group Co_0",
 -202: "Internal Error in operation of group 2^{1+24}.Co_1",
 -203: "Element of 2^{1+24}.Co_1 is not in subroup 2^{1+24}",
 -204: "Bad tag for atom in group 2^{1+24}.Co_1", 

}



def error_string(code):
    """Map an error code to a string"""
    try:
        return QSTATE12_ERROR_STRINGS[code]
    except KeyError:
        err = "Internal error %d in processing QState12 instance"
        return err % code

cpdef int32_t chk_qstate12(int32_t code) except -1:
    """Check the return code of a function in module qstate12.c
    
    The function raises ValueError with a suitable error message if
    the return value ``code``  of a function in module ``qstate12.c``
    is negative. If ``code`` is nonnegative then the function
    returns ``code``.
    """
    if code >= 0:
        return code
    raise ValueError(error_string(code))


####################################################################
# Class QStateMatrix 
####################################################################


cdef class QState12(object):
    """TODO: Yet to be documented!!!"""

    cdef uint64_t data_[QSTATE12_MAXROWS] # data storage
    cdef qstate12_type qs     # structure to store a state vector 
    
    #########################################################################
    # Constructors

    def __cinit__(self, *args, **kwds):
        """Initialize a quadratic state tensor to the scalar 0."""
        cl.qstate12_set_mem(&self.qs, &(self.data_[0]), QSTATE12_MAXROWS);
        cl.qstate12_zero(&self.qs, 0)

    def __init__(self, rows, cols = None, data = None, mode = 0):
        """Intitalize a quadratic state matrix
        
        If ``rows`` is an integer then ``cols`` must also be
        an integer and a quadratic state matrix of shape        
        ``(rows, cols)`` is created. The parameter ``data`` may be:

            * ``None`` (default). Then the zero matrix is created.
            
            * A list of integers. Then that list of integers must 
              encode a valid pair ``(A, Q)`` of bit matrices that make
              up the quadratic state matrix. In this case parameter 
              ``mode`` is  evaluated as follows:
              
               * 1: create matrix ``Q`` from lower triangular part
              
               * 2: create matrix ``Q`` from upper triangular part
               
               * Anything else: matrix ``Q`` must be symmetric.
               
              In this case ``self.factor`` is set to  one.
              
            * An integer ``v``. Then one of the values ``rows``
              and ``cols`` must be zero and the unit (row or
              column) vector with index ``v`` is created.
              This corresponds to a complex matrix with one 
              row or one column.  

        If ``rows`` is an instance of class QState12 then a deep
        copy of that instance is create.                
        """
        cdef uint32_t n, n0, n1, nrows, i
        cdef uint64_t v
        if isinstance(rows, Integral) and isinstance(rows, Integral):
            n0, n1 = rows, cols
            n = n0 + n1
            if data is None:
                chk_qstate12(cl.qstate12_zero(&self.qs, n))                
            elif isinstance(data, Iterable):
                nrows = len(data)
                self.qs.reduced = 0
                for i in range(min(nrows, QSTATE12_MAXROWS)):
                    self.data_[i] = data[i]
                chk_qstate12(cl.qstate12_set(
                    &self.qs, n, nrows, &self.data_[0], mode))
            elif isinstance(data, Integral) and rows * cols == 0:
                v = data
                chk_qstate12(cl.qstate12_vector_state(
                    &self.qs, n, v))
            else:
                err = "Bad data type for QStateMatrix"
                raise TypeError(err) 
            chk_qstate12(cl.qstate12_mat_reshape(&self.qs, n0, n1))
        elif isinstance(rows, QState12) and cols is None:
            source_pqs = pqs12(rows)
            chk_qstate12(cl.qstate12_copy(source_pqs, &self.qs))
        else:
            err = "Cannot construct QStateMatrix from given input" 
            raise TypeError(err) 



    #########################################################################
    # Properties and methods for obtaining implementation details

    property data:
        """Get a copy of the data bit matrix of the state as a list"""
        def __get__(self):
            """Return a copy of the data bits of the  state as a numpy array"""
            chk_qstate12(cl.qstate12_reduce(&self.qs))
            chk_qstate12(cl.qstate12_check(&self.qs))
            cdef uint32_t i
            cdef uint32_t imax = self.qs.nrows
            cdef uint64_t *a = &self.data_[0]
            return  [int(a[i]) for i in range(imax)]
            
    property raw_data:
        """Get a copy of the raw data bit matrix of the state as a list
        
        For debugging and testing only!
        """
        def __get__(self):
            cdef uint32_t i
            cdef uint32_t imax = QSTATE12_MAXROWS
            cdef uint64_t *a = &self.data_[0]
            return  [int(a[i]) for i in range(imax)]
        
    property nrows:
        """Get the number of rows of the bit matrix of the state """
        def __get__(self):
            return self.qs.nrows 

    property ncols:
        """Get the number of columns the bit matrix A of the state 
        
        ``self.ncols == n`` means that the state ``qs`` referred by
        ``self`` is a complex ``2**n`` dimensional vector describing
        a stabilizer state of ``n`` qubits.        
        """
        def __get__(self):
            return self.qs.ncols 

    property factor:
        """Get the factor of the state 
        
        The factor is returned as a pair ``(e, phi)`` which denotes
        the factor :math:`2^{e/2} \cdot \exp(\phi \pi i /4)` .           
        """
        def __get__(self):
            cdef int64_t r = self.qs.nrows
            cdef int64_t f = self.qs.factor
            cdef int64_t e = (f & -0x10) >> 4 if r else 0
            cdef int64_t phi = f & 7 if r else 0
            return (e, phi)


                                
    def row(self, uint32_t i):    
        """Return row i of the bit matrix of the state"""
        self.check()
        cdef uint64_t v
        if i >= self.qs.nrows:
            chk_qstate12(-7)
        v = self.data_[i]
        return int(v)
        
    def column(self, uint32_t j):  
        """Return column j of the bit matrix of the state"""
        self.check()
        return int(cl.qstate12_get_column(&self.qs, j))
            
    def mul_av(self, uint64_t v):
        """Return matrix product ``A * transposed(v)``
        
        Here ``v`` is an integer denoting a bit vector and ``A`` is the 
        ``A`` of the bit matrix of the state.
        
        The result is returned as an integer denoting a bit matrix.
        """
        cdef uint64_t w
        chk_qstate12(cl.qstate12_mul_Av(&self.qs, v, &w))
        return int(w)


    #########################################################################
    # Elementary operations and checks
    
    def is_reduced(self):
        return bool(self.qs.reduced)

    property shape:
        """Get shape of the complex matrix represented by the state 
        
        The function returns a pair ``(rows, cols)`` meaning
        that the state corresponds to a complex
        ``2**nrows`` times ``2**ncols`` matrix.        
        """
        def __get__(self):
            cdef uint32_t cols = self.qs.shape1
            cdef uint32_t rows = self.qs.ncols - cols
            return rows, cols  

    def copy(self):
        return self.__class__(self)

    def conjugate(self):
        """Conjugate state in place"""
        chk_qstate12(cl.qstate12_conjugate(&self.qs))
        return self
        
    def mul_scalar(self, int32_t e, int32_t phi=0):
        """Multiply the state by a sclar factor in place
        
        The factor is :math:`2^{e/2} \cdot \exp(\phi \pi i /4)` .           
        """
        chk_qstate12(cl.qstate12_mul_scalar(&self.qs, e, phi & 7))
        return self
        
    def transpose(self):
        """Transpose state matrix im place"""
        chk_qstate12(cl.qstate12_mat_t(&self.qs))
        return self

    def lb_rank(self):
        """Yet to be documented"""
        cdef int32_t res = cl.qstate12_mat_lb_rank(&self.qs)
        if (res >= -1):
            return res
        chk_qstate12(res)
        
    def set_zero(self):
        """Set state matrix to zero"""
        self.qs.factor = self.qs.nrows = self.reduced = 0
        return self

    def check(self):
        """Raise ValueError is state matrix is bad""" 
        cl.qstate12_check(&self.qs)
        return self

    def check_code(self):
        """Check a state
        
        Return ``0`` if state is ok.
        
        If the state is bad then return an error code as a 
        negative integer.
        """ 
        cdef int32_t res = cl.qstate12_check(&self.qs)
        return res 

    def monomial_row_op(self):
        a = np.zeros(self.qs.nrows, dtype = np.uint32)
        cdef uint32_t[:] a_view = a 
        cdef int32_t res
        res = cl.qstate12_monomial_matrix_row_op(&self.qs, &a_view[0])
        if (res >= 0):
            return a
        chk_qstate12(res)

    #########################################################################
    # Reshaping a state matrix
    
    def reshape(self, int32_t rows, int32_t cols):
        chk_qstate12(cl.qstate12_mat_reshape(&self.qs, rows, cols))
        return self

    #########################################################################
    # Reducing a state matrix to a standard form 

    def echelon(self):
        """Convert state to (not reduced) echelon form

        Only the first  ``n_reduce`` colums are converted to reduced
        echelon form

        The function reutrns a list ``table`` with ``table[j] = i``
        if the leading coefficient of row ``i`` is in columng ``j``.
        We have  ``table[j] = -1``  if no such row exists.
        """
        chk_qstate12(cl.qstate12_echelonize(&self.qs))
        return self

    def reduce(self, forced = False):
        """Reduce a state to a standard form
        
        The representation of a state is not unique.  This function
        changes a sate to a unique standard form where the bit matrix A
        of the internal representation of the state is in reduced
        echelon form, as described in the guide.
        
        If ``forced`` is ``True``, also a previously reduced state
        is recuced one again.
        """
        if forced:
            self.qs.reduced = 0
        chk_qstate12(cl.qstate12_reduce(&self.qs))
        return self 

    def reduce_matrix(self):
        """Yet to be documented"""
        cdef uint8_t row_table[QSTATE12_MAXCOLS+1]
        chk_qstate12(cl.qstate12_reduce_matrix(&self.qs, &row_table[0]))
        cdef uint32_t i = self.nrows + self.ncols
        cdef uint32_t j
        return [ row_table[j] for j in range(i) ]

 
                

    #########################################################################
    # Permuting the qubits of the state

    def rot_bits(self, int32_t rot, uint32_t nrot, uint32_t n0 = 0):  
        """Wrapper for the corresponding C function"""
        chk_qstate12(cl.qstate12_rot_bits(&self.qs, rot, nrot, n0))
        return self
        

    def xch_bits(self, uint32_t sh, uint64_t mask):
        """Wrapper for the corresponding C function"""
        chk_qstate12(cl.qstate12_xch_bits(&self.qs, sh, mask))
        return self

    #########################################################################
    # Applying qubit gates

    def gate_not(self, uint64_t v):
        """Wrapper for the corresponding C function"""
        chk_qstate12(cl.qstate12_gate_not(&self.qs, v))
        return self

    def gate_ctrl_not(self, uint64_t vc, uint64_t v):
        """Wrapper for the corresponding C function"""
        chk_qstate12(cl.qstate12_gate_ctrl_not(&self.qs, vc, v))
        return self

    def gate_phi(self, uint64_t v, uint32_t phi):
        """Wrapper for the corresponding C function"""
        chk_qstate12(cl.qstate12_gate_phi(&self.qs, v, phi))
        return self

    def gate_ctrl_phi(self, uint64_t v1, uint64_t v2):
        """Wrapper for the corresponding C function"""
        chk_qstate12(cl.qstate12_gate_ctrl_phi(&self.qs, v1, v2))
        return self

    def gate_h(self,  uint64_t v):
        """Wrapper for the corresponding C function"""
        chk_qstate12(cl.qstate12_gate_h(&self.qs, v))
        return self
        
    #########################################################################
    # Extending and restricting a state

    def extend_zero(self, uint32_t j, uint32_t nqb):
        """Wrapper for the corresponding C function"""
        chk_qstate12(cl.qstate12_extend_zero(&self.qs, j, nqb))
        return self

    def extend(self, uint32_t j, uint32_t nqb):
        """Wrapper for the corresponding C function"""
        chk_qstate12(cl.qstate12_extend(&self.qs, j, nqb))
        return self
        
    def restrict_zero(self, uint32_t j, uint32_t nqb):
        """Wrapper for the corresponding C function"""
        chk_qstate12(cl.qstate12_restrict_zero(&self.qs, j, nqb))
        return self

    def restrict(self, uint32_t j, uint32_t nqb):
        """Wrapper for the corresponding C function"""
        chk_qstate12(cl.qstate12_restrict(&self.qs, j, nqb))
        return self

    def sumup(self, uint32_t j, uint32_t nqb):
        """Wrapper for the corresponding C function"""
        chk_qstate12(cl.qstate12_sum_cols(&self.qs, j, nqb))
        return self
        
    #########################################################################
    # Conversion of a state matrix to a complex vector
        
    def complex(self):
        """Convert the state to a complex vector"""
        cdef uint32_t n0, n1
        n0, n1 = self.shape
        a = np.empty(2 << self.ncols, dtype = np.double)
        cdef double[:] a_view = a
        chk_qstate12(cl.qstate12_complex(&self.qs, &a_view[0]))
        c = a[0::2] + 1j * a[1::2]
        del a
        return c.reshape((1 << n0, 1 << n1))

    #########################################################################
    # Comparing state matrices

    def __eq__(self, QState12 other):
        """Return True iff two states are equal"""
        cdef p_qstate12_type other_pqs
        if isinstance(other, QState12):
            other_pqs = pqs12(other)
            return chk_qstate12(cl.qstate12_equal(&self.qs, other_pqs))
        err = "Cannot compare QStateMatrix with non QStateMatrix"
        raise TypeError(err)
        
        
    
    #########################################################################
    # Matrix multiplication and inversion
       
    def __imatmul__(self, QState12 other):
        """Perform a matrix multiplication.
        
        We replace ``self`` by the matrix product ``self @ other``.
        ``other`` is not changed.        
        """
        cdef int32_t res
        cdef p_qstate12_type pqs2
        if isinstance(other, QState12):
            pqs2 = pqs12(other)
            res = cl.qstate12_matmul(&self.qs, pqs2, &self.qs)
            if res >= 0:
                return self
            err = "Multiplying QStateMatrix objects of shape %s and %s"
            print("\n" + err % (self.shape, other.shape))
            chk_qstate12(res)
        err = "Bad type for matrix multiplication with QStateMatrix"
        raise TypeError(err)

    def __matmul__(self, QState12 other):
        """Return matrix product  ``self @ other``.
        
        """
        return self.copy().__imatmul__(other)

    def _mat_inv(self):
        """Return inverse matrix"""
        qs1 = self.copy()
        cdef p_qstate12_type pqs1 = pqs12(qs1)
        chk_qstate12(cl.qstate12_mat_inv(pqs1))
        return qs1

    #########################################################################
    # Scalar multiplication 


    def  __imul__(self, other):
        cdef p_qstate12_type pqs
        cdef uint32_t nqb, shape1
        if isinstance(other, Complex):
            qstate12_mul_scalar(self, other)
            return self
        elif isinstance(other, QState12):
            if self.shape == other.shape:
                pqs = pqs12(other)
                nqb = self.ncols
                shape1 = self.qs.shape1
                chk_qstate12(cl.qstate12_product(&self.qs, pqs, nqb, 0))
                self.qs.shape1 = shape1
                return self
            else:
                err = "QStateMatrix instances must have same shape"
                raise ValueError(err)
        else:
            err = "Bad type for multiplication with QStateMatrix"
            raise TypeError(err)


    def  __mul__(self, other):
        return self.copy().__imul__(other)
        
    __rmul__ = __mul__    

    def __itruediv__(self, other):
        if isinstance(other, Complex):
            return self.__imul__(1.0/other)
        else:
            err = "Bad type for multiplication with QStateMatrix"
            raise TypeError(err)

    def __truediv__(self, other):
        return self.copy().__itruediv__(other)
            
    def __neg__(self):
        return self.copy().mul_scalar(0, 4)    

    def __pos__(self):
        return self    

    
    #########################################################################
    # A rather general multiplication function

    def qstate12_product(self, QState12 other, uint32_t nqb, uint32_t nc):
        """Wrapper for the corresponding C function
    
        We return a certain product ``qs3`` of the states ``qs1`` 
        ad ``qs2`` referred by ``self`` and ``other``.
        
        Let ``n1 = qs1.ncols,  n2 = qs2.ncols``.
        Put ``qs1a =  qs1.extend(n1, n2-nqb)``,
        ``qs2a =  qs2.extend(nqb, n1-nqb)``. Then ``qs1a`` and
        ``qs2a`` are complex functions on ``(nn1 + nn2 - nqb)``
        qubits. Let ``qs3a`` be the complex function which is the 
        product of the functions  ``qs1a`` and ``qs2a``. Then we 
        have ``qs3 = qs3a.sum_cols(0, nc)``. So we have
        ``qs3.ncols = nn1 + nn2 - nqb - nc``.
        The result ``qs3`` is computed in the state referred by 
        ``qs1``. It is reduced. 
 
        E.g. ``qs1.product(qs2, nc, nc)`` is the tensor 
        contraction over the first ``nc`` qubits of ``qs1`` and 
        ``qs2``. ``In case qs1.ncols = qs2.ncols = n``, the 
        function ``qs1.product(qs2, 0, n)`` returns the 
        product of ``qs1`` and ``qs2`` (considered as functions); 
        and ``qs1.product(qs2, n, n)`` returns the scalar
        product of ``qs1`` and ``qs2`` (considered as vectors).
 
        ``qs1.product(qs2, 0, n)`` corresponds to the 
        function ``(qs1 (*) qs2)_n``  where ``(*)`` denotes the 
        ``'\odot'`` operation defined in section 
        *Products and tensor products of quadratic mappings*
        of the guide.   
        """
        qs1 = self.copy()
        cdef p_qstate12_type pqs1 = pqs12(qs1)
        cdef p_qstate12_type pqs2 = pqs12(other)
        chk_qstate12(cl.qstate12_product(pqs1, pqs2, nqb, nc))
        return qs1

    #########################################################################
    # Auxiliary function for Matrix multiplication (for testing)

    def qstate12_prep_mul(self, QState12 other, uint32_t nqb):
        """Wrapper for the corresponding C function
    
        To be used for tests only. Function ``qstate12_prep_mul``
        is an auxliary function for function ``qstate12_product``.    
        """
        qs1 = self.copy()
        qs2 = other.copy()
        cdef p_qstate12_type pqs1 = pqs12(qs1)
        cdef p_qstate12_type pqs2 = pqs12(qs2)
        row_pos = chk_qstate12(cl.qstate12_prep_mul(pqs1, pqs2, nqb))
        return row_pos, qs1, qs2
        

    #########################################################################
    # Obtain complex entries of a state
    
        
    def entries(self, indices):
        """Return complex entries of a state.
        
        Here ``indices`` must be a one-dimensional numpy array of 
        indices of ``dtype = np.uint32``. Unused bits of the 
        indices are ignored.
        
        The function returns a complex numpy array ``c`` of the same 
        shape as the input so that ``c[i]`` is the complex value
        of entry ``indices[i]`` of the state.
              
        This function does not change the state. It uses method 
        ``echelonize`` to bring the representation of the
        state to echelon form.
        """
        cdef uint32_t[:] ind_view = indices
        cdef unsigned int n = len(indices)
        a = np.empty(2 * n, dtype = np.double, order = 'C')
        cdef double[:] a_view = a
        if n:
            chk_qstate12(cl.qstate12_entries(
                &self.qs, n, &ind_view[0], &a_view[0]))
        c = a[0::2] + 1j * a[1::2]
        del a
        return c

    #########################################################################
    # Some matrix functions
    
    def pauli_vector(self):
        """TODO: yet to be documented!!!"""
        cdef uint64_t v;
        chk_qstate12(cl.qstate12_pauli_vector(&self.qs, &v))
        return int(v)
        
    def pauli_conjugate(self, v, arg = True):   
        v = np.array(v, dtype = np.uint64, copy=True)
        shape = v.shape
        assert len(shape) <= 1
        v = np.ravel(v)
        if len(v) == 0:
            return []
        cdef uint64_t[:] v_view = v
        chk_qstate12(cl.qstate12_pauli_conjugate(
            &self.qs, len(v), &v_view[0], arg))
        if len(shape):
            return [int(x) for x in v]
        else:
            return int(v[0])
 
####################################################################
# Auxiliary functions  
####################################################################

 
cdef p_qstate12_type pqs12(QState12 state):
    """Return pointer to the structure of a ``state`` 
    
    Here a ``state`` must be of type ``QState12``.
    """
    return &state.qs
  
 
####################################################################
# Convert a complex number to a factor for class QStateMatrix 
####################################################################
    

MIN_ABS = 2.0**-1024
EPS = 1.0e-8 
PI =  3.141592653589793238462643383279
LOG2 = 0.693147180559945309417232121458

       
def qstate12_mul_scalar(QState12 qs, x):
    """Multiply quadratic state matrix with scalar factor
    
    Here ``qs`` is an instance of class ``QState12``
    representing a quadratic state matrix and ``x`` is
    a complex number.
    
    This multiplication succeeds if ``x = 0`` or 
    ``x == 2**(0.5*e) * z**phi``, with ``z = (1+1j)/sqrt(2)``
    and integers ``e, phi``.
    
    Raise ValueError multiplication fails.. We accept a relative 
    error of about ``1.0e-8`` for ``x``.
    
    So it is safe to write e.g. ``2**-0.5 * (-1+1j) * m`` 
    if ``m`` is an  instance of class ``QState12``.    
    """
    cdef double r, phi, e, e_r, phi8_r
    cdef int32_t e_i, phi8_i
    r, phi = cmath.polar(x)
    if r <= MIN_ABS:
        return qs.set_zero()
    e = 2.0 * log(r) / LOG2
    phi8 = 4.0 * phi / PI
    e_r, phi8_r = round(e), round(phi8)
    e_i, phi8_i = <int32_t>e_r, <int32_t>phi8_r
    if (max(fabs(e - e_r), fabs(phi8 - phi8_r)) < EPS):
        return qs.mul_scalar(e_i, phi8_i)
    err = "Cannot convert number to factor for QStateMatrix"
    raise ValueError(err)    




 
####################################################################
# Wrappers for C functions contructing state matrices
####################################################################

def qstate12_unit_matrix(QState12 qs, uint32_t n):
    """Change state qs to a  2**n times 2**n unit matrix"""
    cdef p_qstate12_type m_pqs = pqs12(qs)
    chk_qstate12(cl.qstate12_unit_matrix(m_pqs, n))
    return qs


def qstate12_pauli_matrix(QState12 qs, uint32_t n, uint64_t v):
    """Change state qs to a  2**n times 2**n Pauli matrix"""
    cdef p_qstate12_type m_pqs = pqs12(qs)
    chk_qstate12(cl.qstate12_pauli_matrix(m_pqs, n, v))
    return qs


def qstate12_column_monomial_matrix(QState12 qs, uint32_t nqb, a):
    """Change state ``qs`` to a monomial matrix

    Set the state ``qs`` to a real monomial 
    ``2**nqb`` times ``2**nqb`` transformation matrix ``T`` which 
    is monomial in  the sense that each column contains exactly 
    one nonzero entry ``1`` or ``-1``. So left multiplication  
    with ``T`` maps maps unit vectors to (possibly negated) unit  
    vectors. It transforms a  state ``<v|`` of ``nqb`` input bits 
    to a state ``T * <v|`` of ``nqb`` output qubits. Here
    ``|v>`` is considered as a column vector.

    ``pa`` refers to an array a of integers ``a[i]`` of length 
    ``nqb + 1``. Each integer ``a[i]`` is interpreted as a bit 
    field  via its binary  representation. So ``a[i,j]`` means 
    ``(a[i] >> j) & 1``. ``a[i, j1:j2]`` means the bit field  
    ``a[i,j1],...,a[i,j2-1]``.
   
    For any bit vector ``v`` of length ``nqb`` let ``|v>`` be the 
    unit vector with index ``v``. For any bit vector ``v`` of 
    length ``nqb + 1`` let ``|v>`` be the (possibly negated) unit 
    vector ``(-1)**v[nqb] * |v[0:nqb]>``.  ``|v1 ^ v2>`` and 
    ``|1 << v1>`` are defined via the corrresponding operators 
    ``<<`` and ``^`` in C.
  
    Then ``T``  maps
  
      * ``|0>``     to  ``|a[0, 0:nqb+1]>``
  
      * ``|1 << i>`` to  ``|a[0, 0:nqb+1] ^ a[i+1, 0:nqb+1]>``
  
    ``T`` maps unit vectors to (possibly negated) unit vectors, 
    so ``T(v)`` is well defined by ``|T(v)> = T(|v>)`` for a bit 
    field ``v`` of length ``nqb + 1``. We have
 
      * ``|T(v1 ^ v2)> = (-1)**Q(v1,v2) * |T(v1) ^ T(v2) ^ T(0)>``,
 
    for bit fields ``v1, v2`` of length ``nqb + 1`` and an 
    alternating bilinear form ``Q`` depending on the first ``nqb`` 
    bits of ``v1`` and ``v2`` only. Thus ``T`` is  defined by the 
    above equation and ``Q``. The bilinear form ``Q`` is defined by:
     
      * ``Q(v1, v2) = Q(v2, v1),  Q(v1, v1) = 0``,  and
 
      * ``Q(1 << i, 1 << j) =  m[i + 1, j + nqb + 1]``,  
        for ``0 <= j < i < nqb``.
    """
    cdef p_qstate12_type m_pqs = pqs12(qs)
    cdef uint64_t aa[QSTATE12_MAXROWS+1]
    if nqb >= QSTATE12_MAXROWS:
        return chk_qstate12(-4) 
    cdef int i
    for i in range(nqb + 1):
        aa[i] = a[i]
    chk_qstate12(cl.qstate12_monomial_column_matrix(m_pqs, nqb, &aa[0])) 
    return qs


def qstate12_row_monomial_matrix(QState12 qs, uint32_t nqb, a):
    """Change state ``qs`` to a monomial matrix

    Equivalent to 

    ``qstate12_column_monomial_matrix(qs, nqb, a)``
    ``qs = qs.T``
    """
    cdef p_qstate12_type m_pqs = pqs12(qs)
    cdef uint64_t aa[QSTATE12_MAXROWS+1]
    if nqb >= QSTATE12_MAXROWS:
        return chk_qstate12(-4) 
    cdef int i
    for i in range(nqb + 1):
        aa[i] = a[i]
    chk_qstate12(cl.qstate12_monomial_row_matrix(m_pqs, nqb, &aa[0])) 
    return qs

####################################################################
# Wrappers for other exported C functions  
####################################################################



def qstate12_pauli_vector_mul(uint32_t nqb, uint64_t v1, uint64_t v2):
    if nqb > 31:
        raise ValueError("Bad Pauli group vector")
    return int(cl.qstate12_pauli_vector_mul(nqb, v1, v2))

def qstate12_pauli_vector_exp(uint32_t nqb, uint64_t v, uint32_t e):
    if nqb > 31:
        raise ValueError("Bad Pauli group vector")
    return int(cl.qstate12_pauli_vector_exp(nqb, v, e))


####################################################################
# Support for the subgroup 2^{1+24}.Co_1 of the monster  
####################################################################


@cython.wraparound(False)
@cython.boundscheck(False)
def xsp2co1_chain_short_3(QState12 qstate, src, dest):
    cdef uint64_t[:] src_view = src
    cdef uint64_t[:] dest_view = dest
    cdef p_qstate12_type pqs = pqs12(qstate)
    cdef uint32_t length = len(src)
    assert len(dest) >= length
    if length > 0:
        chk_qstate12(cl.xsp2co1_chain_short_3(pqs, length, &src_view[0],
            &dest_view[0]))

@cython.wraparound(False)
@cython.boundscheck(False)
def xsp2co1_elem_to_qs_i(elem):
    assert len(elem) >= 26 
    cdef uint64_t e[26]
    cdef uint32_t i
    for i in range(26): e[i] = elem[i]
    result = QState12(12, 12)
    cdef qstate12_type qs0
    cdef p_qstate12_type pqs = pqs12(result)
    chk_qstate12(cl.xsp2co1_elem_to_qs_i(e, &qs0))
    chk_qstate12(cl.qstate12_copy(&qs0, pqs))
    return result

@cython.wraparound(False)
@cython.boundscheck(False)
def xsp2co1_elem_to_qs(elem):
    cdef uint64_t[:] elem_view = elem
    result = QState12(12, 12)
    cdef p_qstate12_type pqs = pqs12(result)
    chk_qstate12(cl.xsp2co1_elem_to_qs(&elem_view[0], pqs))
    return result


@cython.wraparound(False)
@cython.boundscheck(False)
def xsp2co1_qs_to_elem_i(QState12 qstate, uint64_t x1):
    cdef p_qstate12_type pqs = pqs12(qstate)
    result = np.zeros(26, dtype = np.uint64)
    cdef uint64_t[:] result_view = result
    chk_qstate12(cl.xsp2co1_qs_to_elem_i(pqs, x1, &result_view[0]))
    return result


####################################################################
# Wrappers for other exported bit matrix functions  
####################################################################


@cython.wraparound(False)
@cython.boundscheck(False)
def bitmatrix64_t(m, uint32_t ncols):
    m1 = np.array(m, dtype = np.uint64)
    m2 = np.zeros(ncols, dtype = np.uint64)
    cdef uint64_t[:] m1_view = m1
    cdef uint64_t[:] m2_view = m2
    if len(m1) and ncols:
        chk_qstate12(cl.bitmatrix64_t(
            &m1_view[0], len(m1), ncols, &m2_view[0]))
    return m2        


@cython.wraparound(False)
@cython.boundscheck(False)
def bitmatrix64_cap_h(m1, m2, uint32_t j0, uint32_t n):
    cdef uint64_t[:] m1_view = m1
    cdef uint64_t[:] m2_view = m2
    cdef uint32_t i1 = len(m1)
    cdef uint32_t i2 = len(m2)
    cdef uint32_t rows1 = cl.bitmatrix64_echelon_h(&m1_view[0], i1, j0, n)
    cdef uint32_t rows2 = cl.bitmatrix64_echelon_h(&m2_view[0], i2, j0, n)
    cdef uint32_t cap = cl.bitmatrix64_cap_h(&m1_view[0], &m2_view[0], 
         i1, i2, j0, n)
    chk_qstate12(cap)
    return rows1 - cap, rows2 - cap


    
