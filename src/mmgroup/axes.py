"""This is yet a stub and will de documented if future releases

The purpose of this module is to give the user access to
representatives of the orbits of the group :math:`G_{x0}` on
the 2A axes in the Monster. 

A similar function for obtaining the orbits on the 2A axes in
the Baby Monster will also be provided. 

Currently this information is available in the submodules of module ``mmgroup.tests.test_axes``. It is yet in a rather messy state.
"""


from mmgroup.tests.axes.get_baby_sample_axes import get_baby_sample_axes
from mmgroup.tests.axes.get_sample_axes import get_sample_axes
from mmgroup.tests.axes.axis import g_central, g_axis, g_axis_opp
from mmgroup.tests.axes.axis import v_axis, v_axis_opp
from mmgroup.tests.axes.axis import Axis, BabyAxis, set_axis_group
from mmgroup.tests.axes.beautify_axes import beautify_axis
from mmgroup.tests.axes.reduce_axis import reduce_axis_G_x0
from mmgroup.tests.axes.reduce_baby_axis import reduce_baby_axis_G_x0

