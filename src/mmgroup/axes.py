"""Orbits of 2A axes under the action of the group :math:`G_{x0}`

The purpose of this module is to give the user access to
representatives of the orbits of the axes in the Monster
under the action of the group :math:`G_{x0}`. 

A similar function for obtaining the orbits on the 2A axes in
the Baby Monster is also provided. 

.. warning::

  The functions in this module are yet experimental and subject
  to change in the future. We recommend to contact the author
  before using them in research projects.

"""

from mmgroup.tests.axes.axis import Axis, BabyAxis, set_axis_group

try:
    from mmgroup.tests.axes.random_axis import RandomAxis
    from mmgroup.tests.axes.random_axis import RandomBabyAxis
    from mmgroup.tests.axes.random_axis import rand_mm_element
    RANDOM_AXES_SUPPORTED = True
except:
    RANDOM_AXES_SUPPORTED = False
    rand_mm_element = None
    raise


try:
    from mmgroup.tests.axes.griess import product_axis_2A
    from mmgroup.tests.axes.griess import Griess, Griess_scalar
except:
    pass
    #raise



