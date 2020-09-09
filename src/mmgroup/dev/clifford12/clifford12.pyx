# cython: language_level=3
"""Simulation of the  Clifford group of up to 12 qubits.

"""

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

from collections.abc import Iterable
from numbers import Integral

import numpy as np
from libc.stdint cimport uint64_t, uint32_t, int32_t, uint8_t, int8_t

include "../pxd_files/clifford12.pxi"
cimport clifford12 as cl
ctypedef qstate12_type *p_qstate12_type
 

# Messages for error codes returned by the functions in module qstate12.c
QSTATE12_ERROR_STRINGS = {
 -1: "No bit with the requested property found in QState12 instance",
 -2: "QState12 instance is inconsistent",
 -3: "Qubit index error at QState12 instance",
 -4: "State is too large for a QState12 instance",
 -5: "Internal buffer overflow in a QState12 instance",
 -6: "Bit matrix part Q of QState12 instance is not symmetric",
 -7: "Bad row index for bit matrix in a QState12 instance",
 -8: "Internal parameter error in a QState12 instance",
 -9: "Overflow in scalar factor in a QState12 instance",
 -11: "QStateMatrix is not invertible",
 -12: "QStateMatrix is not in the Pauli group",
}



def error_string(code):
    """Map an error code to a string"""
    try:
        return QSTATE12_ERROR_STRINGS[code]
    except KeyError:
        err = "Internal error %d in processing QState12 instance"
        return err % code

cdef int32_t chk_qstate12(int32_t code) except -1:
    """Check the return code of a function in module qstate12.c
    
    The function raises ValueError with a suitable error message if
    the return value ``code``  of a function in module ``qstate12.c``
    is negative. If ``code`` is nonnegative then the function
    returns ``code``.
    """
    if code >= 0:
        return code
    raise ValueError(error_string(code))

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

    def __init__(self, source = 0, data = None, mode = 0):
        """Intitalize a state with ``ncols`` qubits

        If ``source`` is given, it may be another state or an integer
        denoting the number of qubits of the new state.

        If ``source`` is an integer then ``self.factor`` is set 
        to ``1`` and  parameter ``data`` may be:

            * ``None`` (default). Then the zero state is created.

            * An integer ``v``. Then the state is set to ``|v>``.
              
            * A list of integers. Then that list of integers must 
              encode a valid pair ``(A, Q)`` of bit matrices that make
              up the state. In this case parameter ``mode`` is
              evaluated as follows:
              
               * 1: create matrix ``Q`` from lower triangular part
              
               * 2: create matrix ``Q`` from upper triangular part
               
               * Anything else: matrix ``Q`` must be symmetric.
        """
        cdef uint64_t ncols, nrows, i
        cdef p_qstate12_type source_pqs
        if isinstance(source, Integral):
            ncols = source
            if data is None:
                chk_qstate12(cl.qstate12_zero(&self.qs, ncols))
            elif isinstance(data, Integral):
                chk_qstate12(cl.qstate12_vector_state(&self.qs, ncols, data))
            elif isinstance(data, Iterable):
                nrows = len(data)
                for i in range(min(nrows, QSTATE12_MAXROWS)):
                    self.data_[i] = data[i]
                chk_qstate12(cl.qstate12_set(
                    &self.qs, ncols, nrows, &self.data_[0], mode))
            else:
                err = "Bad type of data for creating instance of class  QState12"
                raise TypeError(err)
        elif isinstance(source, QState12):
            source_pqs = pqs12(source)
            chk_qstate12(cl.qstate12_copy(source_pqs, &self.qs))
        else:
            err = "Bad type for creating instance of class  QState12"
            raise TypeError(err)

    #########################################################################
    # Properties and methods for obtaining implementation details

    property data:
        """Get a copy of the data bit matrix of the state as a list"""
        def __get__(self):
            """Return a copy of the data bits of the state as a numpy array"""
            chk_qstate12(cl.qstate12_check(&self.qs))
            cdef int imax = self.qs.nrows
            a = np.empty(imax, dtype=np.uint64)
            cdef uint64_t[:] a_view = a
            cdef int32_t i
            for i in range(imax):
                a_view[i] = self.data_[i]
            # Dealing with numpy arrays of type np.uint64 is a nuisance.
            # So we convert is to a list of ints
            return list(map(int, a)) 
            
    property raw_data:
        """Get a copy of the raw data bit matrix of the state as a list
        
        For debugging and testing only!
        """
        def __get__(self):
            return [int(x) for x in self.data_]
        
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
        
    def set_zero(self):
        """Multiply the state by a sclar factor in place
        
        The factor is :math:`2^{e/2} \cdot \exp(\phi \pi i /4)` .           
        """
        self.qs.factor = self.qs.nrows = 0
        return self

    def check(self):
        """Raise ValueError is state is bad""" 
        cdef int32_t res = cl.qstate12_check(&self.qs)
        return chk_qstate12(res)

    def check_code(self):
        """Check a state
        
        Return ``0`` if state is ok.
        
        If the state is bad then return an error code as a 
        negative integer.
        """ 
        cdef int32_t res = cl.qstate12_check(&self.qs)
        return res 


    #########################################################################
    # Reducing a state to a standard form and checking euality of states

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

    def reduce(self):
        """Reduce a state to a standard form
        
        The representation of a state is not unique.  This function
        changes a sate to a unique standard form where the bit matrix A
        of the internal representation of the state is in reduced
        echelon form, as described in the guide.
        """
        chk_qstate12(cl.qstate12_reduce(&self.qs))
        return self 
 
    def __eq__(self, other):
        """Return True iff two states are equal"""
        cdef uint64_t data1[QSTATE12_MAXROWS] 
        cdef uint64_t data2[QSTATE12_MAXROWS] 
        cdef p_qstate12_type other_pqs = pqs12(other)
        cdef qstate12_type qs1
        cdef qstate12_type qs2
        cl.qstate12_set_mem(&qs1, &(data1[0]), QSTATE12_MAXROWS)
        cl.qstate12_set_mem(&qs2, &(data2[0]), QSTATE12_MAXROWS)
        chk_qstate12(cl.qstate12_copy(&self.qs, &qs1))
        chk_qstate12(cl.qstate12_copy(other_pqs, &qs2))
        chk_qstate12(cl.qstate12_reduce(&qs1))
        chk_qstate12(cl.qstate12_reduce(&qs2))
        return chk_qstate12(cl.qstate12_equal(&qs1, &qs2))

    def row_table(self):
        """Return row table
        
        Yet to be documented!!!
        """
        cdef uint8_t table[QSTATE12_MAXCOLS] 
        chk_qstate12(cl.qstate12_echelonize(&self.qs))
        chk_qstate12(cl.qstate12_echelonize(&self.qs))
        
        
                

    #########################################################################
    # Permuting the qubits of the state

    def rot_bits(self, int32_t rot, uint32_t nrot, uint32_t n0 = 0):  
        """Rotate the qubit arguments the state q
        
        For ``n0 <= i < n0 + nrot``  we map bit ``i`` to bit  
        ``n0 + (i + rot) % nrot==. E.g. ``nrot = 3, rot = 1, n0 = 0``
        means bits are mapped as ``0->1, 1->2, 2->0``.
                
        Let ``n = qs.ncols``. Then the function changes 
        ``qs to qs'`` with
        ``qs'(...,x[n0-1],y[0],...,y[nrot-1],x[n0+nrot],...) =``
        ``qs(...,x[n0-1],z[0],...,z[nrot-1],z[n0+nrot],...,),``
        where ``z[j] = y[j - rot (mod 3)]``.
        """
        chk_qstate12(cl.qstate12_rot_bits(&self.qs, rot, nrot, n0))
        return self
        

    def xch_bits(self, uint32_t sh, uint64_t mask):
        """Exchange qubit arguments the state
        
        We exchange argument bit ``j`` with argument bit ``j + sh``
        of the state if bit ``j`` of ``mask`` is set. If bit ``j``
        of ``mask`` is set then bit ``j + sh`` of ``mask`` must not 
        be set. No ``mask`` bit at  position >= ``self.ncols - sh``
        may be set.
        
        E.g.  ``qs.qstate12_xch_bits(1, 0x11)`` changes the 
        state ``qs``  to a statre ``qs'``  with
        ``qs'(x0,x1,x2,x3,x4,x5,x6,...)``  =
        ``qs(x1,x0,x2,x3,x5,x4,x6,...)``.
        """
        chk_qstate12(cl.qstate12_xch_bits(&self.qs, sh, mask))
        return self
        
    #########################################################################
    # Extending and restricting a state

    def extend_zero(self, uint32_t j, uint32_t nqb):
        """Insert ``nqb`` zero qubits at position ``j``. 
        
        Let ``n = self.ncol`` so that the state ``qs`` referred by 
        ``self`` depends on ``n`` qubits. We change ``qs`` to the
        following state ``qs'`` depending on ``n + ncols`` qubits:
        
        ``qs'(x[0],...,x[j-1],y[0],...,y[nqb-1],x[j],...,x[n0-1])`` 
        is equal to ``qs(x[0],...,x[j-1],x[j],...,[n-1])`` if
        ``y[0] = ... = y[nqb-1] = 0`` and equal to zero otherwise.
        So we increment ``self.ncol`` by ``nqb``.
        
        If the input is reduced then the result is also reduced.
        """
        chk_qstate12(cl.qstate12_extend_zero(&self.qs, j, nqb))
        return self

    def extend(self, uint32_t j, uint32_t nqb):
        """Insert ``nqb`` qubits at position ``j``. . 
        
        Let ``n = self.>ncol`` so that the state ``qs`` referred by 
        ``self`` depends on ``n`` qubits. We change ``qs`` to the
        following state ``qs'`` depending on ``n + ncols`` qubits:
        
        ``qs'(x[0],...,x[j-1],y[0],...,y[nqb-1],x[j],...,x[n0-1])`` 
          = ``qs(x[0],...,x[j-1],x[j],...,[n-1])`` .
        So we increment ``self.>ncol`` by ``nqb``.
        
        If the input is reduced then the result is also reduced.
        """
        chk_qstate12(cl.qstate12_extend(&self.qs, j, nqb))
        return self
        
    def sumup(self, uint32_t j, uint32_t nqb):
        """Sum up ``nqb`` qubits  starting at position ``j``. 
        
        Let ``n =self.ncols`` so that the state ``qs`` referred by 
        ``self`` depends on ``n`` qubits. We change ``qs`` to the
        following state ``qs'`` depending on ``n - ncols`` qubits:
        
        ``qs'(x[0],...,x[j-1],x[j+nqb],...,x[n-1])`` =
           ``sum_{{x[j],...,x[j+nqb-1]}}  qs1(x[0],...,x[nn1-1])`` .
        So we decrement ``self.ncol`` by ``nqb``.
        
        The output is not reduced.
        """
        chk_qstate12(cl.qstate12_sum_cols(&self.qs, j, nqb))
        return self

    def restrict_zero(self, uint32_t j, uint32_t nqb):
        """Restrict ``nqb`` qubits starting at postion ``j`` to ``0``.
       
        Let ``n = self.ncols`` so that the state ``qs`` referred by 
        ``self`` depends on ``n`` qubits. We change ``qs`` to the
        following state ``qs'`` depending on ``n`` qubits:
       
        qs'(x[0],...,x[n-1]) is equal to qs(x[0],...,x[n-1]) if 
        x[j] = ... = x[j+nqb-1] = 0 and equal to zero otherwise.
        We do not change ``self.ncols``.
        
        The output is reduced if the input is reduced.
        """
        chk_qstate12(cl.qstate12_restrict_zero(&self.qs, j, nqb))
        return self

    def restrict(self, uint32_t j, uint32_t nqb):
        """This is method ``restrict_zero`` with deletion.
       
        Let ``n = self.ncols`` so that the state ``qs`` referred by 
        ``self`` depends on ``n`` qubits.  We change ``qs`` to the 
        following state ``qs'`` depending  on ``n' = n - nqb`` qubits:

        qs'(x[0],...,x[n'-1]) is equal to 
        qs(x[0],...,x[j-1],0,...,0,x[j],...,x[n'-1]). 
        So we decrement ``self.ncols`` by ``nqb``.
       
        The output is not reduced.
        """
        chk_qstate12(cl.qstate12_restrict(&self.qs, j, nqb))
        return self

     
    #########################################################################
    # Applying qubit gates

    def gate_not(self, uint64_t v):
        """Apply not gates to a state

        Change the state ``qs`` referred by ``selg`` to a state ``qs'`` 
        with ``qs'(x) = qs(x (+) v)``, where ``'(+)'`` is the bitwise 
        xor operation.
        The result is reduced if the input is reduced.
        Computing ``qs.gate_not(1 << j)`` corresponds to
        negating qubit ``j``.
        """
        chk_qstate12(cl.qstate12_gate_not(&self.qs, v))
        return self

    def gate_ctrl_not(self, uint64_t vc, uint64_t v):
        """Apply controlled not gates to a state
        
        Change the state ``qs`` referred by ``self`` to a state ``qs'``
        with ``qs'(x) = qs(x (+) <vc,x> * v)``  where ``'(+)'`` is 
        the bitwise  xor operation, and ``<.,.>`` is the scalar 
        product of bit vectors.
        The result is not reduced.
        Computing ``qs.gate_ctrl_not(1 << jc, 1 << j)``, 
        for ``jc != j``, corresponds to applying a controlled not  
        gate  to qubit ``j``  contolled by qubit ``jc``. 
        This operation is unitary if and only if the scalar
        product of ``j`` and ``jc`` is zero.
        """
        chk_qstate12(cl.qstate12_gate_ctrl_not(&self.qs, vc, v))
        return self

    def gate_phi(self, uint64_t v, uint32_t phi):
        """Apply phase gate to a state
        
        Change the state ``qs`` referred by ``self`` to a state ``qs'``
        with ``qs'(x) = qs(x) * sqrt(-1)**(phi * <v,x>)``, where
        ``<.,.>`` is the scalar product of bit vectors and ``'**'`` 
        denotes exponentiation.
        The result is reduced if the input is reduced.
        Computing ``qs.gate_phi(1 << j, phi)`` 
        corresponds to applying a  phase ``(phi * pi/2)``  gate 
        to qubit ``j``. 
        """
        chk_qstate12(cl.qstate12_gate_phi(&self.qs, v, phi))
        return self

    def gate_ctrl_phi(self, uint64_t v1, uint64_t v2):
        """Apply controlled phase gates to a state
        
        Change the state ``qs`` referred by ``self`` to a state ``qs'``
        with ``qs'(x) = qs(x) * (-1)**(<v1,x>*<v2,x>)``, where 
        ``<.,.>``  is the scalar product of bit vectors and ``'**'``
        denotes exponentiation.
        The result is reduced if the input is reduced.
        Computing ``qs.gate_ctrl_phi(1 << j1, 1 << j2)``
        corresponds to applying a  phase ``pi`` gate to 
        qubit ``j2`` controlled by qubit ``j1``. 
        """
        chk_qstate12(cl.qstate12_gate_ctrl_phi(&self.qs, v1, v2))
        return self

    def gate_h(self,  uint64_t v):
        """Apply Hadamard gates to a state
     
        Apply a Hadamard gate to all qubits ``j`` of the state ``qs``
        (referred by ``self``) with  ``v & (1 << j) == 1``.
        Aplying a Hadamard gate to gate ``j`` changes a state ``qs``
        to a state ``1/sqrt(2) * qs'``, where
        ``qs'(..,x[j-1],x_j,x[j+1],..) = qs(..,x[j-1],0,x[j+1],..)``
          + (-1)**(x_j) * qs(..,x[j-1],1,x[j+1],..)  .
        The result is not reduced.
        """
        chk_qstate12(cl.qstate12_gate_h(&self.qs, v))
        return self
        
    #########################################################################
    # Conversion of a state to a complex vector
    
    
    def complex(self):
        """Convert the state to a complex vector"""
        cdef uint64_t data[QSTATE12_MAXROWS] 
        cdef qstate12_type qs
        cl.qstate12_set_mem(&qs, &(data[0]), QSTATE12_MAXROWS)
        chk_qstate12(cl.qstate12_copy(&self.qs, &qs))
        chk_qstate12(cl.qstate12_reduce(&qs))
        a = np.empty(2 << self.ncols, dtype = np.double)
        cdef double[:] a_view = a
        chk_qstate12(cl.qstate12_complex_test(&qs, &a_view[0]))
        c = a[0::2] + 1j * a[1::2]
        del a
        return c
    
    def complex_unreduced(self):
        """Convert the state to a complex vector without reducing it
       
        To be used for testing only!
        """
        a = np.empty(2 << self.ncols, dtype = np.double)
        cdef double[:] a_view = a
        chk_qstate12(cl.qstate12_complex_test(&self.qs, &a_view[0]))
        c = a[0::2] + 1j * a[1::2]
        del a
        return c

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
    
    def pauli_vector(self, uint32_t n):
        """TODO: yet to be documented!!!"""
        cdef uint64_t v;
        chk_qstate12(cl.qstate12_pauli_vector(&self.qs, n, &v))
        return int(v)
        
    def pauli_conjugate(self, uint32_t n, v):   
        v = np.array(v, dtype = np.uint64, copy=True)
        shape = v.shape
        assert len(shape) <= 1
        np.ravel(v)
        if len(v) == 0:
            return []
        cdef uint64_t[:] v_view = v
        chk_qstate12(cl.qstate12_pauli_conjugate(
            &self.qs, n, len(v), &v_view[0]))
        if len(shape):
            return list(map(int,v))
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
  
def as_qstate12(QState12 state):
    """Return an object as an instance of class QState12"""
    return state
 
 
####################################################################
# Wrappers for exported C functions  
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

def qstate12_mat_t(QState12 qs, uint32_t nqb):
    """Wrapper for the corresponding C function"""
    cdef p_qstate12_type m_pqs = pqs12(qs)
    chk_qstate12(cl.qstate12_mat_t(m_pqs, nqb))


def qstate12_product(QState12 qs1, QState12 qs2, uint32_t nqb, uint32_t nc):
    """Wrapper for the corresponding C function
    
    We compute a certain product ``qs3`` of the states ``qs1`` 
    and store that product in ``qs1``. ``qs2`` is destroyed.
        
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
    cdef p_qstate12_type pqs1 = pqs12(qs1)
    cdef p_qstate12_type pqs2 = pqs12(qs2)
    chk_qstate12(cl.qstate12_product(pqs1, pqs2, nqb, nc))

def qstate12_prep_mul(QState12 qs1, QState12 qs2, uint32_t nqb):
    """Wrapper for the corresponding C function
    
    To be used for tests only. Function ``qstate12_prep_mul``
    is an auxliary function for function ``qstate12_product``.    
    """
    cdef p_qstate12_type pqs1 = pqs12(qs1)
    cdef p_qstate12_type pqs2 = pqs12(qs2)
    return chk_qstate12(cl.qstate12_prep_mul(pqs1, pqs2, nqb))
    
    
def qstate12_matmul(QState12 qs1, QState12 qs2, uint32_t nqb):
    """Perform a matrix multiplication.
        
    Given a state ``qs`` with ``qs->ncols = m + n`` we
    may interpret ``qs`` as an ``2**m`` times ``2**n`` matrix, 
    where the ``2**m`` rows of the matrix correspond to the
    higher ``m`` qubits and the ``2**n`` columns corrspond to 
    the lower ``n`` qubits.
    Then we consider the state ``qs1`` as a ``qs1->ncols - nqb`` 
    times ``nqb`` matrix ``m1``and the state ``qs2`` as an
    ``nqb`` times ``qs2->ncols - nqb`` 
    matrix ``m2``. 
    
    We compute the matrix product ``m3 = m1 @ m2``  as a 
    ``qs1->ncols - nqb`` times ``qs2->ncols - nqb`` matrix
    and we store the state corresponding to matrix ``m3``
    in ``qs1``. `qs2`` is destroyed.
    """
    cdef p_qstate12_type pqs1 = pqs12(qs1)
    cdef p_qstate12_type pqs2 = pqs12(qs2)
    return chk_qstate12(cl.qstate12_matmul(pqs1, pqs2, nqb))

def qstate12_reduce_matrix(QState12 qs, uint32_t nqb):
    """Yet to be documented"""
    cdef uint8_t row_table[QSTATE12_MAXCOLS+1]
    cdef p_qstate12_type pqs = pqs12(qs)
    chk_qstate12(cl.qstate12_reduce_matrix(pqs, nqb, &row_table[0]))
    cdef uint32_t i = qs.nrows + qs.ncols
    cdef uint32_t j
    return [ row_table[j] for j in range(i) ]
    

def qstate12_mat_lb_rank(QState12 qs, uint32_t nqb):
    """Yet to be documented"""
    cdef p_qstate12_type pqs = pqs12(qs)
    cdef int32_t res = cl.qstate12_mat_lb_rank(pqs, nqb)
    if (res < -1):
       chk_qstate12(res)
    return res
    
def qstate12_mat_inv(QState12 qs, uint32_t nqb):
    """Yet to be documented"""
    cdef p_qstate12_type pqs = pqs12(qs)
    chk_qstate12(cl.qstate12_mat_inv(pqs, nqb))
    

def qstate12_bit_matrix_t(m, uint32_t ncols):
    m1 = np.array(m, dtype = np.uint64)
    m2 = np.zeros(ncols, dtype = np.uint64)
    cdef uint64_t[:] m1_view = m1
    cdef uint64_t[:] m2_view = m2
    if len(m1) and ncols:
        chk_qstate12(cl.qstate12_bit_matrix_t(
            &m1_view[0], len(m1), ncols, &m2_view[0]))
    return m2        


def qstate12_pauli_vector_mul(uint32_t nqb, uint64_t v1, uint64_t v2):
    if nqb > 31:
        raise ValueError("Bad Pauli group vector")
    return int(cl.qstate12_pauli_vector_mul(nqb, v1, v2))

def qstate12_pauli_vector_exp(uint32_t nqb, uint64_t v, uint32_t e):
    if nqb > 31:
        raise ValueError("Bad Pauli group vector")
    return int(cl.qstate12_pauli_vector_exp(nqb, v, e))
