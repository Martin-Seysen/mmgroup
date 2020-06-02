
.. include:: definitions.rst

=========================
The mmgroup API reference
=========================


Introduction
============

The *mmgroup* package is a python implementation of Conway's 
construction :cite:`Con85` of the monster group :math:`\mathbb{M}`. 
In mathematics the *monster group* :math:`\mathbb{M}` is the largest 
sporadic finite simple group. The group :math:`\mathbb{M}` has order

   :math:`2^{46} \cdot 3^{20} \cdot 5^9 \cdot 7^6 \cdot 11^2 \cdot 13^3 \cdot 17 \cdot 19 \cdot 23 \cdot 29 \cdot 31 \cdot 41 \cdot 47 \cdot 59 \cdot 71` 
   = :math:`\small 808.017.424.794.512.875.886.459.904.961.710.757.005.754.368.000.000.000` . 

It has a rational representation :math:`\rho` of dimension 
:math:`196884` with coefficients in :math:`\mathbb{Z}[\frac{1}{2}]`,
see :cite:`Con85`. 
So we obtain a modular representation :math:`\rho_p` of 
:math:`\mathbb{M}` by reducing the rational representation 
:math:`\rho` modulo an arbitrary odd number :math:`p`. The current
version of the *mmgroup* package supports the representations 
:math:`\rho_p` of :math:`\mathbb{M}` for 
:math:`p = 3, 7, 15, 31, 127, 255`.

The *mmgroup* package uses highly optimized C programs for 
calculating in the representations of :math:`\rho_p` of
:math:`\mathbb{M}`, which are most efficient if :math:`p+1` 
is a small power of two. Most of these C programs have been
generated automatically by a program code generator which is
also part of this project, but not part of the public 
interface described in this document.

In the description of the *mmgroup* package we use the
notation in :cite:`Seysen20`, where an explicit description
of the representations :math:`\rho(g), g \in G` is given
for a generating subset  :math:`G` of 
:math:`\mathbb{M}`.


Installation
============

.. include:: installation.inc

Basic structures
================

.. _golay-label:

The Golay code and its cocode
-----------------------------

.. automodule:: mmgroup.structures.gcode


.. autoclass:: mmgroup.GCode
   :members:

.. autoclass:: mmgroup.GcVector
   :members:

.. autoclass:: mmgroup.Cocode
   :members:

.. autoclass:: mmgroup.Parity
   :members:



.. _parker-loop-label:


The Parker loop
---------------


.. automodule:: mmgroup.structures.ploop

.. autoclass:: mmgroup.PLoop
   :members:

.. autofunction:: mmgroup.PLoopZ


Octads and suboctads
--------------------

.. automodule:: mmgroup.structures.suboctad

.. autofunction:: mmgroup.Octad

.. autoclass:: mmgroup.SubOctad
   :members:

.. _aut_ploop_label:

Automorphisms of the Parker loop
--------------------------------

.. automodule:: mmgroup.structures.autpl


.. autoclass:: mmgroup.AutPL
   :members: 


.. _basis-golay-label:

The basis of the Golay code and of its cocode
---------------------------------------------


.. include:: golay_basis.inc



.. _mmgroup-label:

The Monster Group
=================


.. automodule:: mmgroup.mm_group


.. autoclass:: mmgroup.MMGroup
   :members: atom, neutral, rand_word, rand_atom, 
             sample, word,   


.. autoclass:: mmgroup.MMGroupWord
   :members:  as_tuples, copy,  is_reduced, reduce, order 


.. _mmrep-label:

The Representation of the Monster Group
=======================================


.. automodule:: mmgroup.mm_space


.. autoclass:: mmgroup.MMSpace
   :members:  zero, rand_uniform, from_sparse, from_bytes,
              tuple_to_index, index_to_tuple, index_to_short

.. autoclass:: mmgroup.MMSpaceVector
   :members:  as_bytes, as_sparse, as_tuples, projection,
              mul_exp


.. autofunction:: mmgroup.characteristics



.. only:: html

   .. rubric:: **References**

.. bibliography:: references.bib



