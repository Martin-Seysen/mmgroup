import numpy as np

from mmgroup import MM0, mat24 
from mmgroup.tests.axes.get_sample_axes import get_sample_axes

from mmgroup.tests.axes.equations import solve_Qx0_equations

#######################################################################

NONZERO_ROWS = []
OCCUR = {}
LEN_OCCUR = 3
MAX_OCCUR = 64


def nonzero_rows_6F():
    global NONZERO_ROWS    
    axis = get_sample_axes()['6F']
    for i in range(759):
        v = mat24.octad_to_vect(i)
        if v & 0xff == 0:
            xrow = (axis)['T', i]
            s_occur = sum(1 << k for k in range(64) if xrow[k])
            if s_occur:        
                NONZERO_ROWS.append(i)
            if len(NONZERO_ROWS) >= LEN_OCCUR:
                break

initialized = False


def analyse_6F_test_elements():
    yield MM0()
    data = set([0, 255])
    for octad in range(759):
        gc = mat24.octad_to_gcode(octad)
        v = mat24.gcode_to_vect(gc) & 0xff
        if v not in data:
             yield MM0('y', gc)
             data |= set([v, 255 - v])


def analyse_cases_6F(verbose = 0):
    dict_occur = {}
    axis = get_sample_axes()['6F']
    for n, g in enumerate(analyse_6F_test_elements()):
        ax = axis * g
        if verbose > 2:
            print("Analyzing test case", n)            
            ax.display_sym(0, text = "A")
            ax.display_sym('B', text = "B")
            ax.display_sym('C', text = "C")
        occur = []
        for i in NONZERO_ROWS:
            xrow = (ax)['T', i]
            s_occur = sum(1 << k for k in range(64) if xrow[k])
            assert s_occur
            occur.append(s_occur)
        assert len(occur) == LEN_OCCUR
        occur = tuple(occur)
        if occur not in dict_occur:
            dict_occur[occur] = str(g**-1)
        if len(dict_occur) >= MAX_OCCUR:
            break
    return dict_occur


def analyse_6F(verbose = 0, reanalyse = False):
    global initialized
    if initialized and not reanalyse:
        return
    if verbose:
        print("Analysing orbit 6F...")
    global OCCUR
    nonzero_rows_6F()
    OCCUR = analyse_cases_6F(verbose)
    assert len(OCCUR) == MAX_OCCUR,  len(OCCUR)
    if verbose:
        print("Orbit 6F analysed, %d cases found" % len(OCCUR))
        print("Rows of tag T to be tested:", NONZERO_ROWS)
    if verbose > 1:
        for n, (occ, g) in enumerate(OCCUR.items()):
            print("OCCUR entry %d, g = %s" % (n, g))
            for i, s in enumerate(occ):
                print("%2d, 0x%016x" % (NONZERO_ROWS[i], s))
    initialized = True


def solve_6F(ax):
    occ = []
    for i in NONZERO_ROWS[:LEN_OCCUR]:
        row = (ax)['T', i]
        occ.append(sum(1 << k for k in range(64) if row[k]))
    occ = tuple(occ)
    if occ in OCCUR:
       return OCCUR[occ]


def solve_special_equation_6F(axis):
    if not initialized:
        analyse_6F(verbose = 0)
    ax = axis.copy()
    g0 = axis.group(solve_6F(ax))
    ax *= g0
    g1 = solve_Qx0_equations('6F', ax)
    return g0 * g1


def test_find_6F(verbose = 0):
    orbit = '6F'
    print("Test axis reduction in orbit 6F")
    ref_axis = get_sample_axes()['6F'] 
    for n in range(300):
        if verbose > 1:
            print("Test case", n) 
        ax = ref_axis * MM0('r', 'G_x0')
        ax =  beautify_axis(ax.copy())
        g = solve_special_equation_6F(ax)
        ax *= g
        if verbose > 2:
            ax.display_sym(0, text = "A")
            ax.display_sym('B', text = "B")
            ax.display_sym('C', text = "C")
        for tag in "ABCTXYZ":
            if (ax[tag] != ref_axis[tag]).any():
                print("Axes differ in tag", tag)
        assert ax.v15 == ref_axis.v15
    print("Test passed")
    

if __name__ == "__main__":
    VERBOSE = 1
    analyse_6F(VERBOSE)
    from mmgroup.tests.axes.beautify_axes import beautify_axis
    from mat24_orbits import initialize_all
    initialize_all()
    test_find_6F(VERBOSE)


  


