

================================
The C interface
================================


Introduction
============

This *document* describes the functionality of the C modules in this
project. For most of these C modules there is also a python extension.




Description of the ``mmgroup.mat24`` extension
==============================================

.. automodule:: mmgroup.dev.mat24.mat24_doc


C interface
-----------

.. doxygenfile:: mat24_functions.c




Generating C code for the ``mmgroup.mat24`` extension
-----------------------------------------------------

In this section we give a brief overview over the modules in
``mmgroup.dev.mat24`` used for generating C code for the 
``mmgroup.mat24`` extension.



Module ``mat24_ref``
....................

.. automodule:: mmgroup.dev.mat24.mat24_ref

Module ``mat24tables``
......................

.. automodule:: mmgroup.dev.mat24.mat24tables

Module ``mat24aux``
...................

.. automodule:: mmgroup.dev.mat24.mat24aux


Module ``mat24heptad``
......................

.. automodule:: mmgroup.dev.mat24.mat24heptad

Module ``make_addition_table``
..............................

.. automodule:: mmgroup.dev.mat24.make_addition_table

Module ``make_mul_transp``
..........................

.. automodule:: mmgroup.dev.mat24.make_mul_transp


Description of the ``mmgroup.mat24_xi`` extension
=================================================

This extension is used for computing tables containing the monomial
part of the operation of the generators :math:`\xi^i, i=1,2`.
These elements operate monomially on the ``98280`` entries of the 
representation :math:`\rho_p` with tags ``B``, ``C``, ``T``, and ``X``. 
One table is generated for each exponent :math:`i`.

Module ``mmgroup.dev.mat24_xi.mat24_xi_ref`` is a pure python 
substitute for this extension; but calculating the tables for the
generators :math:`\xi^i` with that extension takes a rather long time.

More details will be documented in a future version of this project.


Description of the ``mmgroup.mm`` extension
===========================================


Module ``mmgroup.mm`` provides basic functions for the representation
of the monster. 

Module ``mmgroup.mm`` is implemented as an extension with ``Cython``.
The main source file for that extension is ``mm_basics.pyx`` in
directory ``src.mmgroup.dev.mm_basics``. Most functions exported by
this extension are in one-to-one correspondence with certain C
functions exported from several automatically generated .c files. 

Representation of a vector in :math:`\rho_p`
--------------------------------------------

.. automodule:: mmgroup.dev.mm_basics.mm_doc

Basic table-providing classes for module ``mmgroup.mm``
-------------------------------------------------------

.. automodule:: mmgroup.dev.mm_basics.mm_basics

.. autoclass:: mmgroup.dev.mm_basics.mm_basics.MM_Const
   :members: snippet



More details will be documented in a future version of this project.


Description of the ``mmgroup.mm<p>`` extensions
===============================================


Module ``mmgroup.mm<p>`` implements the representation :math:`\rho_p`
for a fixed modulus :math:`p = 2^k-1`. E.g  module ``mmgroup.mm3`` 
implements :math:`\rho_3`. 

Module ``mmgroup.mm<p>`` is implemented as an extension with 
``Cython``. The main source file for the ``mmgroup.mm<p>`` extension
is the automatically generated file ``mm_op<p>.pyx`` in directory 
``src.mmgroup.dev.pxd_files``. Each function exported by such an 
extension corresponds to a certain C function  exported from an 
automatically generated .c file.

For reasons of efficiency a dedicated set of C files is generated for 
each modulus ``p``. Corresponding C files for different moduli are 
generated from the same source file with extension ``.ske``. 


The basic table-providing class for ``mmgroup.mm<p>``
-----------------------------------------------------

.. automodule:: mmgroup.dev.mm_op.mm_op

.. autoclass:: mmgroup.dev.mm_op.mm_op.MM_Op


More details will be documented in a future version of this project.


Description of the ``qstate12`` and ``clifford12`` extensions
=============================================================

C functions dealing with quadratic state vectors
................................................


The C functions in modules ``qstate.c`` and ``qsmatrix.c`` perform 
operations on quadratic state matrices given by triples 
:math:`(e, A, Q)` as described in **API reference** the section 
:ref:`clifford-group-label`. 


Check everything from here on!!!!

A quadratic state vector :math:`v` of type
``qbstate12_type`` with component ``ncols = n`` models a complex 
vector in a vector space  :math:`V` of dimension :math:`2^n`, and 
that the basis of ``V`` is labelled by the elements of the Boolean
vector space :math:`\mathbb{F}_2^n`. In C and python programs
we represent the element :math:`(x_{n-1}, \ldots, x_{0})` of
:math:`\mathbb{F}_2^n` by the integer 
:math:`\sum_{0 \leq i < n} 2^i \cdot x_i`. This leads to a natural
representation of ``v`` as a one-dimensional complex array of
length :math:`2^n`, starting with index ``0``.

A quadratic state matrix is a quadratic shape vector augmented
by an information about its matrix shape. For a quadratic state
matrix ``qs`` of shape ``(n0, n1)`` we put ``qs.ncols = n0 + n1``
and ``qs.ncols = n1``.

The zero state is encoded as a matrix with :math:`m'=0` rows.
We do not update the entries :math:`Q_{i,0}`, so the 
corresponding bits in compoment ``data`` of the structure
are garbage. One may use the C function ``qstate12_check`` to
set these bits to their proper values.


The current implementation requires ``n + m <= 63``.  
This can easily be generalized to larger Clifford 
groups by reserving an array of several integers for each row 
of matrix :math:`M`. Here we also leave the details to the reader.

C functions supporting this module are prefixed with ``qbstate12_``.
Unless otherwise stated, these functions return an ``int32_t``, 
where a nonegative value is interpreted as success, and a negative 
value is intepreted as failure. Depending on the function, a 
nonnegative return value may e.g. mean an index for a matrix
:math:`A`, :math:`M`, or :math:`Q`.

Typical names for parameters of functions in this module are:

   ================== ================================================
   ``pqs, pqs1, ...`` Pointer to structure of type ``qbstate12_type``
   ``nqb``            Number of qubits, i.e. of columns of matrix 
                      :math:`A`.
   ``nrows``          Number of rows of matrix :math:`A`, :math:`M`, 
                      and :math:`Q`.
   ``i, i1, ...``     Index of a row of matrix :math:`A`, :math:`M`,  
                      or and :math:`Q`, starting with 0.
   ``j, j1, ...``     Index of a column of matrix :math:`A`, with a 
                      column of :math:`A`, corrsesponding to a qubit, 
                      starting with ``j = 0``.
                      If appropriate, an index  ``j >= ncols`` refers 
                      to column ``(j - ncols)`` of matrix  :math:`Q`.
   ``pv, pv1,...``    Pointer to a row or column vector of matrix 
                      :math:`A`, :math:`M`, or  :math:`Q`.
   ================== ================================================
   
   .. doxygenfile:: clifford12.h

   .. doxygenfile:: qstate12.c




