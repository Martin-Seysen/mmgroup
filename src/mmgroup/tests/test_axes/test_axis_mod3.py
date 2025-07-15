import sys
from random import randint
import numpy as np
import pytest

if __name__ == "__main__":
    sys.path.append("../../../")


from mmgroup import MM0, MMV, AutPL
V = MMV(255)
PI_XCHG_23 = MM0(AutPL(0, {1:1, 2:3, 3:2}, False))
MODE_GROUPS = {0: 'M', 1: 'B', 2: 'B & 2E_6'}


def import_all():
    global Axis, BabyAxis, mm_profile_mod3_permute24
    from mmgroup.axes import Axis, BabyAxis
    from mmgroup.mm_reduce import mm_profile_mod3_permute24

def rand_Nxyz(mode = 0):
    r"""Return random *even* element of the group :math:`N_{x0}`

    The function returns a random even element :math:`g` of the group
    :math:`N_{x0}` as an instance of class `MM0`. Element :math:`g`
    also fixes a certain set :math:`\{x_\delta \mid \delta \in S\}` of
    2A involutions pointwise, where :math:`S` is a set of Golay cocode
    words depending of parameter ``mode``. Legal values for ``mode``
    are:

    === ====================================
     0  :math:`S = \{ \}` 
     1  :math:`S = \{ [2,3] \}`
     2  :math:`S = \{ [1,2], [2,3] \}`
    === ====================================
    """
    g = MM0('r', 'N_x0_e & %s' % MODE_GROUPS[mode])
    if mode == 2:
        # Permutation pi in M_24 corresponding to g fixes 1; but it may
        # either fix or exchange 2 and 3. Change g so that pi fixes 2
        # and 3, if necessary.
        pi = g.as_M24_permutation()
        assert pi[1] == 1
        if pi[2] == 3:
            g *= PI_XCHG_23
            pi = g.as_M24_permutation()
        assert pi[1:4] == [1,2,3]
    return g



def display_matrix(m, text = "", skip_zero = False):
    print(text)
    for j1 in range(24):
        print("%1d:" % (j1 % 8), end = " ") 
        for j2 in range(24):           
            if skip_zero and m[j1, j2] == 0:
                print("  .", end = " ")
            else:
                print("%03x" % m[j1, j2],  end = " ")
            if j2 & 7 == 7: print("", end = " ")
        print("")
    print("")

def display_hash_A(h, show_sort = False):
    _, hash, a = h
    display_matrix(a)
    print("hash = %016x" % hash)



def ax_hash_equal(hash1, hash2, verbose = 0):
    mat1, mat2 = hash1[2], hash2[2]
    equ_mat = (mat1 == mat2).all()
    equ = equ_mat and hash1[1] == hash2[1]
    if equ and not verbose:
        return True
    if verbose or not equ:
        MAX_ERR = 10
        if  hash1[1] != hash2[1]:
            print("Axes have different hash values for N_x0_e mod 2")
        print("Hash matrices are %sequal" % ("" if equ_mat else "un"))
        if (verbose > 1 or not equ_mat):     
            display_hash_A(hash1, 0) 
            display_hash_A(hash2, 0)
        if not equ_mat: 
            display_matrix(mat1 ^ mat2, "Matrix difference", 1)
    elif verbose:
        print("Axes have same hash value for N_x0_e mod 2")
    return equ     


def permute_matrix24(b, g):
    b1 = np.zeros(576, dtype = np.uint16)
    try:
        g = np.array(g.as_M24_permutation(), dtype = np.uint8) 
    except:
        g = np.array(g, dtype = np.uint8)
    pi = np.array(g, dtype = np.uint8)
    b = np.array(b,  dtype = np.uint16)
    assert b.shape in [(576,), (24,24)]
    assert mm_profile_mod3_permute24(b.ravel(), pi, b1) == 0
    return b1.reshape((24,24))


def do_test_mm_profile_mod3_permute24(verbose = 0):
    import_all()
    for i in range(2):
        g = MM0('r', 'AutPL')
        v = V('R')
        m1, m2 = (v * g)['A'],  permute_matrix24(v['A'].ravel(), g)
        if verbose:
            print("Testing function mm_profile_mod3_permute24")
            if verbose > 1:
                display_matrix(v['A'], "A")
                display_matrix(m1, "m1")
                display_matrix(m2, "m2")
        assert (m1 == m2).all()  



def one_test_axis_profile(mode = 0, test_S3 = True, verbose = 1):
    ax = Axis('r')
    ax_hash = ax.profile_Nxyz(mode = mode)
    for j in range(2):
        if verbose:
            print("\nTest N_xyz", j, ", mode =", mode)
        g = rand_Nxyz(mode)
        ax1 = ax * g
        ax1_hash = ax1.profile_Nxyz(mode = mode)
        #ax1_hash[0][13,17] += 1 # cause a bug
        if verbose:
            print("g =", g)
            print(g.as_M24_permutation())
        equ = ax_hash_equal(ax1_hash, ax_hash, verbose)
        assert equ
        permuted = permute_matrix24(ax_hash[0], g)
        ok = (ax1_hash[0] == permuted).all()
        if verbose or not ok:
            err = "Hash value is not permutation invariant"
            print(err)
            display_matrix(ax1_hash[2], "Hash")
            display_matrix(ax1_hash[0], "Matrix")
            print("Permutation\n", g.as_M24_permutation())
            display_matrix(permuted, "Permuted")
            display_matrix(ax1_hash[0] ^ permuted, "Diff", True)
            if not ok:
                raise ValueError(err)

    if test_S3:
        for e in range(3):
            for f in range(2):
                if verbose:
                    print("\nTest S_3, mode =", mode, ", S_3=", e, f)
                h1 = ax1.profile_Nxyz((e,f), mode = mode)
                ax_m = ax * MM0([('t',e), ('d', 0x800 & -f)])
                h2 = ax_m.profile_Nxyz(mode = mode)
                assert ax_hash_equal(h1, h2, verbose)






@pytest.mark.axes
def test_axis_profile(n_axes = 5, verbose = 0):
    import_all()
    do_test_mm_profile_mod3_permute24(verbose)
    for i in range(n_axes):
        for mode in [0, 1, 2]:
            if verbose:
                print("\nTest", i+1, ", mode =", mode)
            one_test_axis_profile(mode, i < 3, verbose = verbose)
           

# Axis().g_axis * Axis(MM(AX_4B)).g_axis  is of class 4B
AX_4B = "y_324h*x_17c8h*d_122h*p_69815652*l_1*p_2027520*l_1*p_593649*l_1*t_1*l_2*p_2956800*l_1*p_21507363*l_2*t_2*l_2*p_1985280*l_1*p_85374263*l_1*p_11616000"

@pytest.mark.axes
def test_axis_product_class(n_axes = 5, verbose = 0):
    import_all()
    SQ_DICT = {
    "2A":"1A", "2B":"1A", "4A":"2B", "4B":"2A", "4C":"2B", "6A":"3A",
    "6C":"3A", "6F":"3C", "8B":"4A", "10A":"5A", "10B":"5A", "12C":"6A"
    }
    x = MM0('x', 0x1000)
    return
    for cl, ax in {**Axis.representatives(), **BabyAxis.representatives()}.items():
        if cl[-1:].isdigit(): cl = cl[:-1]
        sp = ax.product_class(ax * x)
        assert sp == SQ_DICT[cl], (cl, sp, SQ_DICT)
        ax1 = ax.copy()
        sp1 = ax.product_class(ax * x, sparse = True)
        assert sp1 == sp

    assert Axis().product_class(Axis(MM(AX_4B))) == "4B"



 
@pytest.mark.slow
@pytest.mark.bench
@pytest.mark.axes
def test_bechmark_axis_profile():
    import_all()
    import time
    NTESTS = 10000
    from mmgroup.mm_reduce import mm_profile_mod3_load
    from mmgroup.mm_reduce import mm_profile_mod3_hash
    a = np.zeros((64,72), dtype = np.uint64)
    b = np.zeros(2*24*24, dtype = np.uint16) 
    for i in range(len(a)):
        mm_profile_mod3_load(15, Axis('r').v15.data, a[i], 0)
    t_start = time.process_time()
    for i in range(NTESTS):
        mm_profile_mod3_hash(a[i & 63], b, 0)
    t = time.process_time() - t_start
    print("\nRunime for watermarking N_x0 orbit of axis: %.4f ms" %
        (t * 1000.0 / NTESTS))
 


