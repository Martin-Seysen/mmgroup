

======================================
The C interface of the mmgroup project
======================================


Introduction
============

This *document* describes the functionality of the C modules in this
project. For most of these C modules there is also a python extension.


.. highlight:: none


Description of the ``mmgroup.mat24`` extension
==============================================

.. automodule:: mmgroup.dev.mat24.mat24_doc


C interface
-----------


.. doxygenfile:: mat24_functions.h

.. doxygenfile:: mat24_functions.c

.. doxygenfile:: mat24_random.c


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


.. _mmgroup-generators-label:

Description of the ``mmgroup.generators`` extension
===================================================

.. automodule:: mmgroup.dev.generators.generators_doc


.. _header-mmgroup-generators-label:


Header file mmgroup_generators.h
--------------------------------

.. doxygenfile::  mmgroup_generators.h
 

Implementing the group :math:`N_{0}`
------------------------------------


.. automodule:: mmgroup.dev.generators.mm_group_n_doc

C interface for file mm_group_n.c
..................................

.. doxygenfile:: mm_group_n.c


.. _op-leech-label:

Operation of the group :math:`G_{x0}` on the Leech lattice
----------------------------------------------------------


.. automodule:: mmgroup.dev.generators.gen_leech_doc


C interface for file gen_leech.c
.................................

.. doxygenfile:: gen_leech.c


C interface for file gen_leech_type.c
......................................

.. doxygenfile:: gen_leech_type.c


C interface for file gen_leech3.c
.................................

.. doxygenfile:: gen_leech3.c


C interface for file gen_leech_reduce.c
........................................

.. doxygenfile:: gen_leech_reduce.c


C interface for file gen_leech_reduce_n.c
.........................................

.. doxygenfile:: gen_leech_reduce_n.c


Support for the generator  :math:`\xi` of the monster group
-----------------------------------------------------------   

.. automodule:: mmgroup.dev.generators.gen_xi_ref

C interface for file gen_xi_functions.c
........................................

.. doxygenfile:: gen_xi_functions.c


A random generator
------------------

.. automodule:: mmgroup.dev.generators.gen_random_doc



C interface for file gen_random.c
..................................

.. doxygenfile:: gen_random.c





Description of the ``qstate12`` and ``clifford12`` extensions
=============================================================

C functions dealing with quadratic state vectors
------------------------------------------------


The C functions in modules ``qstate.c`` and ``qsmatrix.c`` perform 
operations on quadratic state matrices given by triples 
:math:`(e, A, Q)` as described in **API reference** the section 
:ref:`clifford-group-label`. 



A quadratic state vector :math:`v` of type
``qbstate12_type`` with component ``ncols = n`` models a complex 
vector in a vector space  :math:`V` of dimension :math:`2^n`, 
or an  :math:`2^{n-k} \times 2^k` matrix.

The basis of vector space :math:`V` is labelled by the elements of 
the  Boolean vector space :math:`\mathbb{F}_2^n`. In C and python 
programs we represent the element :math:`(x_{n-1}, \ldots, x_{0})` 
of :math:`\mathbb{F}_2^n` by the integer 
:math:`\sum_{0 \leq i < n} 2^i \cdot x_i`. This leads to a natural
representation of a vector :math:`v \in V` as a one-dimensional 
complex array of length :math:`2^n`, starting with index 
:math:`0`.

A quadratic state matrix is a quadratic shape vector augmented
by an information about its matrix shape. For more details we 
refer to the description of struct ``qbstate12_type`` in file
``clifford12.h``.

  
Header file ``clifford12.h``
----------------------------
   .. doxygenfile:: clifford12.h

Module ``bitmatrix64.c``
------------------------

   .. doxygenfile:: bitmatrix64.c

Module ``uint_sort.c``
------------------------

   .. doxygenfile:: uint_sort.c


Module ``qstate12.c``
---------------------
     
   .. doxygenfile:: qstate12.c

Module ``qstate12io.c``
------------------------
     
   .. doxygenfile:: qstate12io.c

Module ``qmatrix12.c``
----------------------
     
   .. doxygenfile:: qmatrix12.c


.. _c-functions-G-x0-label:

C functions dealing with the subgroup :math:`G_{x0}` of the monster
-------------------------------------------------------------------

.. automodule:: mmgroup.dev.clifford12.xsp2co1_doc



Module ``xsp2co1.c``
----------------------

   .. doxygenfile:: xsp2co1.c


Module ``xsp2co1_elem.c``
--------------------------

   .. doxygenfile:: xsp2co1_elem.c


Module ``leech3matrix.c``
--------------------------

   .. doxygenfile:: leech3matrix.c


Module ``involutions.c``
--------------------------

   .. doxygenfile:: involutions.c


Module ``xsp2co1_traces.c``
---------------------------

   .. doxygenfile:: xsp2co1_traces.c


Module ``xsp2co1_map.c``
--------------------------

   .. doxygenfile:: xsp2co1_map.c




Description of the ``mmgroup.mm`` extension
===========================================


Module ``mmgroup.mm`` provides basic functions for the representations
of the monster modulo small numbers ``p``.  For each supported
modulus ``p`` there is also a specific module ``mmgroup.mm<p>`` 
containing highly optimized C programs dealing with the
representation modulo ``p``.

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



Header file mm_basics.h
--------------------------------

.. doxygenfile::  mm_basics.h
 

C interface for file mm_aux.c
---------------------------------------

.. doxygenfile:: mm_aux.c


C interface for file mm_tables.c
---------------------------------------

.. doxygenfile:: mm_tables.c

C interface for file mm_random.c
---------------------------------------

Module ``mm_random.c`` is deprecated an no longer in use.


Description of the ``mmgroup.mm<p>`` extensions
===============================================


Module ``mmgroup.mm<p>`` implements the 196884-dimensional rational
representation :math:`\rho_p` of the monster group modulo a fixed 
odd modulus :math:`p = 2^k-1`. E.g  module ``mmgroup.mm3`` 
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


Naming conventions
------------------

For each supported characteristic ``p`` there is a cython extension
``mm<p>`` in module ``mmgroup``. So e.g. for ``p = 15`` there is the
extension ``mm15``.

For any such extension the code generator generates several C files 
from each of the source files with extension ``.ske`` stored in 
subdirectory ``mmgroup.dec.mm_op``. So e.g. from the source file
``mm_op_pi.ske`` we generate the file ``mm3_op_pi.c`` for ``p=3``,
``mm7_op_pi.c`` for ``p=7``, etc. The code generator provides
different instances of class ``mmgroup.dev.mm_op.MM_Op`` for different 
characteristics ``p``. Any of these classes provides certain macro
variables and functions for the corresponding characteristic ``p``.
Here the most important variable is ``P``, which is set to the
characteristic ``p``. So e.g. for ``p = 31`` the code generator
changes every occurrence of the string ``%{P}`` in a ``.ske`` file
to the string ``31`` in the generated C file. Of course, there are 
also macros for more sophisticated replacements.

Every exported C function in one of the modules ``mm<p>`` has a name
starting with ``mm_op%{P}_`` in the ``.ske`` file containing its
implementation. The code generation process also generates a 
``.pyx`` file with name ``mm_op<p>.pyx`` for each supported 
characteristic ``p``. In the ``.pyx`` file a C function with
prefix ``mm_op%{P}_`` is wrapped into a Cython function with
prefix ``op_``.

So e.g. in case ``p = 15``, for the C function with name 
``mm_op%{P}_pi`` in file ``mm_op_pi.ske`` we generate a C function
with name ``mm_op15_pi`` in file ``mm15_op_pi.c``. We also
generate the file ``mm_op15.pyx`` containing a function with name
``op_pi`` that wraps the C function ``mm15_op_pi.c``. It is an
(arguable) advantage that cython functions doing the same job
for different characteristics ``p`` have the same name.

A C function and its cython wrapper have the same parameters.
Here all parameters are integers or  pointers to integers.
A python function usually passes a ``numpy`` array (of
appropriate type) as an argument to a cython function, where 
the corresponding C function expects a pointer to an integer. 
The return type of such a C function is either ``void`` or
some integer type.

In the following subsections of this section we will document 
the C functions for characteristic ``p = 15`` only.

Special functions in characteristic ``p = 15`` and ``p = 3``
............................................................

In principle, the Cython extensions ``mmgroup.mm<p>`` support
the same set of functions for all characteristics ``p``.

In a later phase of the project it has turned out that support 
for dealing with certain special vectors in the representation 
:math:`\rho_p` may speed up the operation in the monster group 
by a factor of about 1000! These special vectors are called 2A  
axes in :cite:`Con85`.

The implementation of the additional functions for supporting 2A 
axes is a bit involved; and it turns out that it suffices to 
deal with the case ``p = 15``. So we add some extra functions 
to the Cython extension ``mmgroup.mm15``. Some of these
function are also added to the extension ``mmgroup.mm3``.

This comprises modules ``mm15_op_rank_A.c``, ``mm15_op_eval_A.c``, 
and ``mm15_op_eval_X.c``. A module ``mm3_op_rank_A.c`` (for ``p = 3``)
is also present.



The basic table-providing class for ``mmgroup.mm<p>``
-----------------------------------------------------

.. automodule:: mmgroup.dev.mm_op.mm_op

.. autoclass:: mmgroup.dev.mm_op.mm_op.MM_Op



Header file mm_op15.h
--------------------------------

.. doxygenfile::  mm_op15.h


C interface for file mm15_op_misc.c
---------------------------------------

.. doxygenfile:: mm15_op_misc.c


C interface for file mm15_op_pi.c
---------------------------------------

.. doxygenfile:: mm15_op_pi.c


C interface for file mm15_op_xy.c
---------------------------------------

.. doxygenfile:: mm15_op_xy.c



C interface for file mm15_op_t.c
---------------------------------------

.. doxygenfile:: mm15_op_t.c

C interface for file mm15_op_xi.c
---------------------------------------

.. doxygenfile:: mm15_op_xi.c


C interface for file mm15_op_word.c
---------------------------------------

.. doxygenfile:: mm15_op_word.c



C interface for file mm15_op_rank_A.c
---------------------------------------

.. doxygenfile:: mm15_op_rank_A.c



C interface for file mm15_op_eval_A.c
---------------------------------------

.. doxygenfile:: mm15_op_eval_A.c


C interface for file mm15_op_eval_X.c
---------------------------------------

.. doxygenfile:: mm15_op_eval_X.c



Description of the ``mm_reduce`` extension
==========================================

The functions in this module implement the fast reduction of an
element of the monster group described in :cite:`Seysen22`. There
we define a triple of vector :math:`(v_1, v^+, v^-)` in the
representation :math:`\rho_{15}` such that an element 
:math:`g` of the monster van be recognized from the triple  
:math:`(v_1 \cdot g, v^+ \cdot g, v^- \cdot g)`. A precomputed
vector :math:`v_1` is stored in file ``mm_order_vector.c``.
Module ``mm_order.c`` contains functions for computing the
order of an element of the monster. Module ``mm_reduce.c``
contains the function ``mm_reduce_M`` that implements the
fast reduction algorithm in the monster group.
 


Generating an **order vector**
------------------------------

Python module ``mmgroup.dev.mm_reduce.find_order_vector``
.........................................................

.. automodule:: mmgroup.dev.mm_reduce.find_order_vector


Header file mm_reduce.h
--------------------------------------

.. doxygenfile:: mm_reduce.h



C interface for file mm_order_vector.c
--------------------------------------

.. doxygenfile:: mm_order_vector.c


C interface for file mm_order.c
-------------------------------

.. doxygenfile:: mm_order.c



C interface for file mm_reduce.c
--------------------------------

.. doxygenfile:: mm_reduce.c


C interface for file mm_shorten.c
---------------------------------

.. doxygenfile:: mm_shorten.c


C interface for file mm_vector_v1_mod3.c
-----------------------------------------

.. doxygenfile:: mm_vector_v1_mod3.c




Shared libraries and dependencies between Cython extensions
===========================================================

For each Cython extension the relevant C functions are 
implemented in a shared library (or a DLL in Windows). Names
of the corresponding DLLs in windows are given in the following 
table:


.. table:: Shared libraries implementing Cython extensions
    :widths: 35 65

    ====================== =======================================
    Cython extension       Implemented in shared library
    ====================== =======================================
    ``mm_reduce``          ``mmgroup_reduce``
    ``mm<p>``              ``mmgroup_mm_op<p>``
    ``mm``                 ``mmgroup_mm_basics``
    ``clifford12``         ``mmgroup_clifford12``
    ``generators, mat24``  ``mmgroup_mat24``
    ====================== =======================================

Here the string ``<p>`` is to be replaced by one of the integer 
literals ``3``, ``7``, ``15``, ``31``, or ``255`` to obtain the 
name of the shared library (or of the Cython extension) dealing 
with the representation :math:`\rho_p` of the monster group.
In Linux, the name of each shared library is prefixed by the
string ``lib``. So the Windows DLL  ``mmgroup_mm_op15.dll``
corresponds to the Linux shared library ``libmmgroup_mm_op15.so.``



The Cython extension ``mm_<p>`` depends on other Cython extensions
as shown in the following paragraph:


.. code-block:: none
    
    mm_reduce
       + mm<p>
          + mm
             + clifford12
                + generators
                    + mat24
 




