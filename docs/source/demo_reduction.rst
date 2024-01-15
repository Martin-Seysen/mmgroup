.. include:: definitions.rst


.. _demonstration-label:


===============================================
Demonstration code for the reduction algorithm
===============================================

This chapter is yet under construction!


This chapter has been written for a mathmatically-inclinded reader who
wants to read an implementation demonstrating the reduction algorithm
for the Monster group in :cite:`Seysen22`. 
In the section we present a python implementation of that reduction
algorithm. This implementation has been optimized for readability, 
and not for speed; but it can still be executed and tested.

Our goal is that the reader may acquire a satisfactory demonstration
of the new reduction algorithm by browsing through the python
modules in the ``mmgoup.demo`` package only.

Demonstrating every tiny detail of the reduction algorithm in Python
leads to a bloated software package that is hardly readable or 
executable. So we have to make a compromise, which parts of the
reduction algorithm we want to demonstrate. 

Following the *Pacific island model* in :cite:`Wilson13`, we assume that
computations in the subgroup :math:`G_{x0}` of the Monster (of 
structure :math:`2^{1+24}.\mbox{Co}_1`) are easy enough, so that
demonstration code is not needed for them. The relevant computations
in :math:`G_{x0}` are described in detail in the appendices of
:cite:`Seysen22`. 

We also assume that functions for the operation of the Monster on 
its 196883-dimensional representation :math:`\rho_{15}` (with
coefficients taken modulo 15) are available. This operation
is decscibed in detail in :cite:`Seysen20`. In our demonstration
we use multiplication for that operation.

We also assume that functions performing linear algebra with matrices
over the Leech lattice (modulo 2 and 3) are available.

Module ``mmgroup.demo.redcue_sub`` contains python functions
implementing the required functionality mentioned above. Most
functions in this module are just Python wrappers for the
corresponding C functions in the ``mmgroup`` package.


.. warning::
    Functions and classes presented in Chapter :ref:`demonstration-label`
    are for demonstration only. They should not be used in a
    life application!



Data structures
----------------

.. autoclass:: mmgroup.demo.Mm
   :members: count_triality_elements

.. autoclass:: mmgroup.demo.MmV15


.. autoclass:: mmgroup.demo.Leech2
   :members: type



.. _demo_reduce_sub_label:


Python module reduce_sub.py
---------------------------


.. automodule:: mmgroup.demo.reduce_sub


.. autofunction:: mmgroup.demo.reduce_sub.mat15_norm

.. autofunction:: mmgroup.demo.reduce_sub.mat15_apply

.. autofunction:: mmgroup.demo.reduce_sub.mat15_rank_3

.. autofunction:: mmgroup.demo.reduce_sub.vect15_S

.. autofunction:: mmgroup.demo.reduce_sub.leech2_span

.. autofunction:: mmgroup.demo.reduce_sub.leech2_rad

.. autofunction:: mmgroup.demo.reduce_sub.map_type4_to_Omega

.. autofunction:: mmgroup.demo.reduce_sub.map_type2_to_standard

.. autofunction:: mmgroup.demo.reduce_sub.map_feasible_type2_to_standard

.. autofunction:: mmgroup.demo.reduce_sub.find_triality_element_for_axis

.. autofunction:: mmgroup.demo.reduce_sub.find_in_Nx0

