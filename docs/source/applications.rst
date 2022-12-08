

=======================================
Applications of the **mmgroup** package
=======================================


Introduction
============

Some applications using the *mmgroup* package can be downloaded from
the github repository

https://github.com/Martin-Seysen/mmgroup .

In that repository each application is stored in a subdirectory of 
directory ``applications``.

The *mmgroup* package must be installed before running any of these
applications, as described in the **API reference** of this project.


Representing the Monster as a Hurwitz group
===========================================

.. automodule:: applications.Hurwitz.readme


Mapping the Coxeter group :math:`Y_{555}` into the Bimonster
=============================================================

.. automodule:: applications.Y555.readme


Implementing the Bimonster and the homomorphism from :math:`Y_{555}` to it
---------------------------------------------------------------------------

module Y555.inc_P3
...................

.. automodule:: applications.Y555.inc_p3


.. autoclass:: applications.Y555.inc_p3.P3_node
   :members:  y_name


.. autoclass:: applications.Y555.inc_p3.AutP3
   :members:  order, map, point_map, line_map


module Y555.p3_to_mm
......................


.. automodule:: applications.Y555.p3_to_mm

.. autofunction:: applications.Y555.p3_to_mm.Norton_generators


module Y555.bimm
......................

.. automodule:: applications.Y555.bimm


.. autoclass:: applications.Y555.bimm.BiMM
   :members:  order

.. autofunction:: applications.Y555.bimm.P3_BiMM


.. autofunction:: applications.Y555.bimm.AutP3_BiMM




