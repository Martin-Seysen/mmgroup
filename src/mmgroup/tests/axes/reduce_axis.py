import numpy as np

from mmgroup import MM0, XLeech2, leech2_orbits_raw, mat24, MM
from mmgroup.tests.axes.axis import Axis
from mmgroup.tests.axes.get_sample_axes import get_sample_axes
from mmgroup.tests.axes.beautify_axes import beautify_axis
from mmgroup.tests.axes.equations import solve_Qx0_equations
from mmgroup.tests.axes.case_6F import solve_special_equation_6F
from mmgroup.tests.axes.case_6F import analyse_6F



##################################################################

initialized = False

def reduce_axis_G_x0(ax, check = True, verbose = False):
    global initialized
    orbit = ax.axis_type()
    img_ax = beautify_axis(ax.copy())
    if orbit == '6F':
        if not initialized:
            analyse_6F()
            initialized = True
        img_ax *= solve_special_equation_6F(img_ax)
    elif orbit not in ['2A', '2B']:
        img_ax *= solve_Qx0_equations(orbit, img_ax)
    if check:
        ref = get_sample_axes()[orbit]
        eq = ref.v15 == img_ax.v15
        if (not eq):
            for e in "ABCTXYZ":
                if not (ref[e] == img_ax[e]).all():
                    print("axes differ in tag", e)
            raise ValueError("Reduction of axis in G_x0 failed")
    return  img_ax.g1



