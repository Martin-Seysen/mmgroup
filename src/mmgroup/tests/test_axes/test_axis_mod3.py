import sys
from random import randint
import numpy as np
import pytest

if __name__ == "__main__":
    sys.path.append("../../../")


from mmgroup import MM0, MMV
V = MMV(255)

def import_all():
    global Axis, mm_profile_mod3_permute24
    from mmgroup.axes import Axis
    from mmgroup.mm_reduce import mm_profile_mod3_permute24

def display_hash_A(h, show_sort = False):
    _, hash, a = h
    for i in range(24):
        print("%1d:" % (i % 8), end = " ") 
        for j in range(24): 
            print("%02x" % int(a[i,j]), end = " ")
            if j & 7 == 7: print("", end = " ")
        if i > 0 and show_sort and list(a[i])  < list(a[i-1]):
            for k in range(24):
                if a[i, k] < a[i-1, k]:
                    break
            print("f%d" % k, end = "")
        print("")
    print("hash = %016x" % hash)


def display_matrix(m, text = "", skip_zero = False):
    print(text)
    for j1 in range(24):
        print("%1d:" % (j1 % 8), end = " ") 
        for j2 in range(24):           
            if skip_zero and m[j1, j2] == 0:
                print(" .", end = " ")
            else:
                print("%02x" % m[j1, j2],  end = " ")
            if j2 & 7 == 7: print("", end = " ")
        print("")
    print("")

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
    b1 = np.zeros(576, dtype = np.uint8)
    try:
        g = np.array(g.as_M24_permutation(), dtype = np.uint8) 
    except:
        g = np.array(g, dtype = np.uint8)
    pi = np.array(g, dtype = np.uint8)
    b = np.array(b,  dtype = np.uint8)
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


MODE_GROUPS = {0: 'M', 1: 'B', 2: 'B & 2E_6'}

def one_test_axis_profile(mode = 0, test_S3 = True, verbose = 1):
    ax = Axis('r')
    ax_hash = ax.profile_Nx0(mode = mode)
    for j in range(2):
        if verbose:
            print("\nTest N_x0_e", j, ", mode =", mode)
        g = MM0('r', 'N_x0_e & %s' % MODE_GROUPS[mode])
        ax1 = ax * g
        ax1_hash = ax1.profile_Nx0(mode = mode)
        #ax1_hash[0][13,17] += 1 # cause a bug
        if verbose:
            print("g =", g)
            print(g.as_M24_permutation())
        equ = ax_hash_equal(ax1_hash, ax_hash, verbose)
        assert equ
        assert (ax1_hash[0] == permute_matrix24(ax_hash[0], g)).all()
    if not test_S3:
        return
    for e in range(3):
        for f in range(2):
            if verbose:
                print("\nTest S_3, mode =", mode, ", S_3=", e, f)
            h1 = ax1.profile_Nx0((e,f), mode = mode)
            ax_m = ax * MM0([('t',e), ('d', 0x800 & -f)])
            h2 = ax_m.profile_Nx0(mode = mode)
            assert ax_hash_equal(h1, h2, verbose)






@pytest.mark.axes
def test_axis_profile(n_axes = 10, verbose = 0):
    import_all()
    do_test_mm_profile_mod3_permute24(verbose)
    for i in range(n_axes):
        for mode in [0, 1]: # mode 2 yet buggy!!!
            if verbose:
                print("\nTest", i+1, ", mode =", mode)
                one_test_axis_profile(mode, i < 3, verbose = verbose)
           
            
 
