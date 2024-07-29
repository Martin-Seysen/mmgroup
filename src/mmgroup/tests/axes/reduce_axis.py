import numpy as np

from mmgroup import MM0, XLeech2, leech2_orbits_raw, mat24, MM
from mmgroup.axes import get_sample_axes, Axis
from mmgroup.axes import g_central, g_axis
from mmgroup.axes import beautify_axis
from mmgroup.tests.axes.equations import solve_Qx0_equations
from mmgroup.tests.axes.case_6F import solve_special_equation_6F
from mmgroup.tests.axes.case_6F import analyse_6F



##################################################################

initialized = False

def reduce_axis_G_x0(ax, check = True, verbose = False):
    global initialized
    if not initialized:
        analyse_6F()
        initialized = True
    orbit = ax.axis_type()
    g2 =  beautify_axis(ax.axis_type(), ax.g).g1    
    img_ax = ax * g2
    if orbit in ['2A', '2B']:
        g3 = ax.group()
    elif orbit == '6F':
        g3 = solve_special_equation_6F(img_ax)
    else:
        g3 = solve_Qx0_equations(orbit, img_ax)

    if check:
        img_ax *= g3;
        ref = get_sample_axes()[orbit]
        eq = ref.v15 == img_ax.v15
        if (not eq):
            for e in "ABCTXYZ":
                if not (ref[e] == img_ax[e]).all():
                    print("axes differ in tag", e)
            raise ValueError("Reduction of axis in G_x0 failed")
    return  g2 * g3



