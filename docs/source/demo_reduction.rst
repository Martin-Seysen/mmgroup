.. include:: definitions.rst


.. _demonstration-label:


===============================================
Demonstration code for the reduction algorithm
===============================================


This chapter has been written for a mathmatically-inclinded reader who
wants to read an implementation demonstrating the reduction algorithm
for the Monster group in :cite:`Seysen22`. 
In this section we present a Python implementation of that reduction
algorithm. This implementation has been optimized for readability, 
and not for speed; but it can still be executed and tested.

Our goal is present a satisfactory demonstration of the new reduction
algorithm in this chapter of the project documentation. Executable
Python code can be found in the ``mmgoup.demo`` package. Function
``reduce_monster_element`` in the ``mmgroup.demo.reduce_monster``
packgage actually reduces an element of the Monster.

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
is decscibed in detail in :cite:`Seysen20`. In the ``mmgroup``
package we use multiplication for that operation.

We also assume that functions performing linear algebra with matrices
over the Leech lattice (modulo 2 and 3) are available.


Section :ref:`demonstration-subfunction-label` describes the
interface of functions implementing the tasks mentioned above.




.. warning::
    Functions and classes presented in Chapter :ref:`demonstration-label`
    are for demonstration only. They should not be used in a
    life application!

This demonstration makes some random choices in the reduction
process. This is appropriate for testing. The life implementation
of the reduction algorithm is deterministic.


Data structures
----------------

This section contains the data structures required for the
demonstration of the reduction algorithm. These data structures
are classes that can be imported from module **mmgroup.demo**.


.. autoclass:: mmgroup.demo.Mm
   :members: count_triality_elements

.. autoclass:: mmgroup.demo.MmV15


.. autoclass:: mmgroup.demo.Leech2
   :members: type



.. _demo_reduce_sub_label:



The reduction algorithm for the Monster group
-------------------------------------------------------

.. automodule:: mmgroup.demo.reduce_monster



Module **mmgroup.demo.reduce_monster** requires the following
external references.

.. code-block:: python

    from mmgroup.demo import Mm, Leech2, MmV15
    from mmgroup.demo.reduce_sub import *
    from mmgroup.demo.reduce_axis import reduce_axis
    from mmgroup.demo.reduce_feasible import reduce_feasible_axis


Here comes the main function **reduce_monster_element** of the
module that reduces an element of the Monster group.

.. literalinclude:: ../../src/mmgroup/demo/reduce_monster.py
   :language: python
   :pyobject: reduce_monster_element


Here is a simple test function for testing the main function 
**reduce_monster_element**.


.. literalinclude:: ../../src/mmgroup/demo/reduce_monster.py
   :language: python
   :pyobject: check_reduce_monster_element


The following function reduces an element of the subgroup
:math:`G_{x0}` of the Monster. The implementation of this
function is discussed in :cite:`Seysen22`, Section 6.

.. literalinclude:: ../../src/mmgroup/demo/reduce_monster.py
   :language: python
   :pyobject: reduce_G_x0





Reducing an axis in the Monster
------------------------------------------


.. automodule:: mmgroup.demo.reduce_axis


Module **mmgroup.demo.reduce_axis** requires the following external
references.



.. code-block:: python

    from random import choice                   # returns a random entry of a list
    from mmgroup.demo import Mm, Leech2, MmV15  # data strucures used
    from mmgroup.demo.reduce_sub import *       # functions used

Here is the implementation of function **reduce_axis**.

.. literalinclude:: ../../src/mmgroup/demo/reduce_axis.py
   :language: python
   :pyobject: reduce_axis

Dictionary **TARGET_ORBITS** encodes the graph shown in Figure 2
in :cite:`Seysen22`, Section 8.3. The vertices of that graph are
orbits of axes.

.. code-block:: python

   TARGET_ORBITS = {
     '2B' : ['2A'],
     '4A' : ['2A'],
     '4B' : ['2B'],
     '4C' : ['2B'],
     '6A' : ['4A'],
     '6C' : ['4A'],
     '6F' : ['4C'],
     '8B' : ['4A'],
     '10A' : ['6A'],
     '10B' : ['4B', '4C'],
     '12C' : ['4B', '6A'],
    }

Here is the implementation of function **axis_orbit**.
This function has been implemented according to the discussion of 
the orbits of axes in :cite:`Seysen22`, Section 8.4.


.. literalinclude:: ../../src/mmgroup/demo/reduce_axis.py
   :language: python
   :pyobject: axis_orbit


Here is the implementation of function **compute_U**. For an axis
**v** the function computes the set **U(v)** of vectors in the 
Leech lattice (mod 2) defined in :cite:`Seysen22`,
Section 8.4. The implementation of the function is essentially a 
rewrite of the discussion of the different cases of axis orbits in
that section.


.. literalinclude:: ../../src/mmgroup/demo/reduce_axis.py
   :language: python
   :pyobject: compute_U



Reducing a feasible axis in the Baby Monster
------------------------------------------------


.. automodule:: mmgroup.demo.reduce_feasible


Module **mmgroup.demo.reduce_feasible** requires the 
following external references.


.. code-block:: python


    from mmgroup.demo import Mm, Leech2, MmV15
    from mmgroup.demo.reduce_axis import axis_orbit
    from mmgroup.demo.reduce_axis import compute_U
    from mmgroup.demo.reduce_axis import TARGET_ORBITS
    from mmgroup.demo.reduce_sub import*

Here is the implementation of function **reduce_feasible_axis**.


.. literalinclude:: ../../src/mmgroup/demo/reduce_feasible.py
   :language: python
   :pyobject: reduce_feasible_axis



.. _demonstration-subfunction-label:


Subfunctions for the reduction algorithm
--------------------------------------------

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







