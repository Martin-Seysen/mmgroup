import sys
from random import randint
import pytest

if __name__ == "__main__":
    sys.path.append("../../../")


from mmgroup import MM0


def import_all():
    global Axis
    from mmgroup.axes import Axis

def display_hash_A(h, show_sort = False):
    _, hash, a = h
    for i in range(24):
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

def ax_hash_equal(hash1, hash2, verbose = 0):
    mat1, mat2 = hash1[2], hash2[2]
    equ_mat = (mat1 == mat2).all()
    equ = equ_mat and hash1[1] == hash2[1]
    if equ and not verbose:
        return True
    if verbose or not equ:
        MAX_ERR = 10
        print("Axes have different hash values for N_x0_e mod 2")
        print("Hash matrices are %sequal" % ("" if equ_mat else "un"))
        if not equ_mat:
            done = 0
            for i in range(24):
                for j in range(24):
                    if mat1[i,j] != mat2[i,j]:
                        print("Matrices differ at entry %d, %d (%02x, %02x)"
                             % (i,j, mat1[i,j], mat[i,j]))
                        done += 1
                        if done > MAX_ERR:
                           break
                if done > MAX_ERR:
                    break     
        display_hash_A(hash1, 1) 
        display_hash_A(hash2, 1)
        if not equ_mat: 
            print("Matrix difference")
            display_hash_A((mat1 ^ mat2, 0))
    elif verbose:
        print("Axes have same hash value for N_x0_e mod 2")
    return equ     

@pytest.mark.axes
def test_cases(n_axes = 3, verbose = 0):
    import_all()
    for i in range(n_axes):
        ax = Axis('r')
        ax_hash = ax.profile_Nx0()
        for j in range(2):
            if verbose:
                 print("\nTest", i, j)
            g = MM0('r', 'N_x0_e')
            ax1 = ax * g
            ax1_hash = ax1.profile_Nx0()
            #ax1_hash[0][13,17] += 1 # cause a bug
            if verbose:
                print(g)
            equ = ax_hash_equal(ax1_hash, ax_hash, verbose)
            assert equ
        for e in range(3):
            for f in range(2):
                h1 = ax1.profile_Nx0((e,f))
                ax_m = ax * MM0([('t',e), ('d', 0x800 & -f)])
                h2 = ax_m.profile_Nx0()
                assert ax_hash_equal(h1, h2, verbose)

 
