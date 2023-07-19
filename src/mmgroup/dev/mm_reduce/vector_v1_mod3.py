import os
import sys
import time
import numpy as np
import re
from collections import OrderedDict
from random import randint, sample



if __name__ == "__main__":
   sys.path.append(os.path.join('..', '..', '..'))



import_pending = True


def import_all():
    global import_pending
    if not import_pending: return

    global  MMV, MMVector, Xsp2_Co1, MMV3
    global  vect_to_cocode
    global  ploop_theta
    global  MM0
    global  tuple_to_sparse

    global  gen_leech3to2_type4
    global  leech2matrix_add_eqn
    global  bitmatrix64_t
    global  leech2matrix_solve_eqn
    global  xsp2co1_elem_read_mod3
    global  uint64_parity
    global  gen_leech2_reduce_type4
    global  mm_aux_mmv_extract_sparse_signs
    global  mm_aux_index_sparse_to_leech2
    global  mm_aux_mmv_extract_sparse_signs
    global  mm_aux_mmv_extract_x_signs
    global  mm_aux_mmv_add_sparse
    global  mm_op_eval_A_rank_mod3 
    global  mm_op_watermark_A
    global  mm_op_watermark_A_perm_num
    global  mm_op_word_tag_A
    global  mm_op_checkzero

    from mmgroup import MMV, MMVector, Xsp2_Co1
    from mmgroup.mat24 import vect_to_cocode
    from mmgroup.mat24 import ploop_theta
    from mmgroup.structures.mm0_group import MM0
    from mmgroup.structures.mm_space_indices import tuple_to_sparse

    from mmgroup.generators import gen_leech3to2_type4
    from mmgroup.clifford12 import leech2matrix_add_eqn
    from mmgroup.clifford12 import bitmatrix64_t
    from mmgroup.clifford12 import leech2matrix_solve_eqn
    from mmgroup.clifford12 import xsp2co1_elem_read_mod3
    from mmgroup.clifford12 import uint64_parity
    from mmgroup.generators import gen_leech2_reduce_type4
    from mmgroup.mm_op import mm_aux_mmv_extract_sparse_signs
    from mmgroup.mm_op import mm_aux_index_sparse_to_leech2
    from mmgroup.mm_op import mm_aux_mmv_extract_sparse_signs
    from mmgroup.mm_op import mm_aux_mmv_extract_x_signs
    from mmgroup.mm_op import mm_aux_mmv_add_sparse
    from mmgroup.mm_op import mm_op_eval_A_rank_mod3
    from mmgroup.mm_op import mm_op_watermark_A
    from mmgroup.mm_op import mm_op_watermark_A_perm_num
    from mmgroup.mm_op import mm_op_word_tag_A
    from mmgroup.mm_op import mm_op_checkzero

    MMV3 = MMV(3)




_DIR = os.path.split(__file__)[0]
PY_FILENAME = os.path.join(_DIR, "v1_mod3_data.py")

#######################################################################
# Find a suitable vector v1 (mod 3)
#######################################################################


def make_cocode_indices():
    indices = []
    solve_yt = np.zeros(11, dtype = np.uint64)
    nrows = 0
    for i in range(1, 23):
        eqn = vect_to_cocode(1 + (1 << i))
        found = leech2matrix_add_eqn(solve_yt, nrows, 11, eqn)
        if (found):
            indices.append(i)
            nrows += 1
            if nrows == 11:
                return indices
    raise ValueError("Generation of cocode indices has failed")
        


def r():
    return randint(1,2)


def A_indices():
    indices = make_cocode_indices()
    eqn_indices = [(2, 'A', 0, i) for i in indices]
    fixed_indices = [(2, 'A', i, i) for i in range(22)
        if not i in indices]
    fixed_indices.append( (1, 'A', 22, 22) )
    assert not 22 in indices and not 23 in indices
    rnd_indices = [(r(), 'A', 0, 0)]
    for n, i in enumerate(indices):
         for j in indices[:n+1]:
             rnd_indices.append( (r(), 'A', i, j) )
    rnd_indices = sample(rnd_indices, 13)
    return eqn_indices,  fixed_indices + rnd_indices


def X_indices():
    indices = make_cocode_indices()
    X_indices = [(2, 'B', 0, i) for i in indices]
    X_indices += [(2, 'C', 0, indices[0])]
    x_i = [0] + [1 << i for i in range(11)]
    X_indices += [(2, 'X', i, 0) for i in x_i]
    X_indices += [(2, 'Z', 0, 2), (1, 'Z', 0, 3)]
    return X_indices



def rand_indices():
    rand_i = sample(list(range(1,2047)), 10)
    return [(r(), 'Y', i, randint(0, 23)) for i in rand_i]



def str_watermark(v):
    if isinstance(v, MMVector):
        v = v.data
    wmark = np.zeros(24, dtype = np.uint32)
    if mm_op_watermark_A(3, v.data, wmark) < 0:
        return "Watermarking of vector mod 3 failed"
    else:
        return "Watermark: "+ ",".join(
            ["%06x" % wmark[i] for i in range(9)])    



def check_v1_mod3(v3):
    wmark = np.zeros(24, dtype = np.uint32)
    rank = mm_op_eval_A_rank_mod3(3, v3.data, 0) >> 48
    if rank != 23:
        return False
    if mm_op_watermark_A(3, v3.data, wmark) < 0:
        return False
    eqn_indices, x_indices = A_indices()[0], X_indices()
    for _, t, i, j in eqn_indices + x_indices[:-1]:
        if v3[t, i, j] != 2:
           return False
    for _, t, i, j in x_indices[-1:]:
        if v3[t, i, j] != 1:
           return False
    return True


def make_v1_mod3(verbose = 0):
    ok = False
    wmark = np.zeros(24, dtype = np.uint32)
    for i in range(10000):
       eqn_indices, more_indices = A_indices()
       x_indices = X_indices()
       all_indices = eqn_indices + x_indices + more_indices
       all_indices += rand_indices()
       v3 = MMV3(all_indices)
       if not check_v1_mod3(v3):
           continue
       if verbose:
           print("Compute vector v1 mod 3; No of trials:", i+1)
           rank = mm_op_eval_A_rank_mod3(3, v3.data, 0) >> 48
           print("Rank of part 'A' is", rank)
           print(str_watermark(v3))
       sparse_indices = np.concatenate(
          [tuple_to_sparse(3, *j) for j in all_indices])
       return sparse_indices, v3
    raise ValueError("Could not generate a suitable vector v1 mod 3")
   


HEADER = """# This file has been created automatically, do not change!

# It contains the data of a vector ``v1`` (modulo 3, in the sense
# :cite:`Seysen22`) in sparse representation.

"""


def str_data(text, data):
    s = "%s = [\n   " % text
    for i, x in enumerate(data):
        s += hex(x) + ","
        s += "\n   " if i % 6 == 5 else " "
    s += "\n]\n"
    return s


def write_v1_mod3(result):
    print("Writing file " + PY_FILENAME)
    f = open(PY_FILENAME, "wt")
    print(HEADER, file = f)
    print(str_data("V1_MOD3", list(result)), file = f)
    f.close()
    



def get_v1_mod3_data(recompute = False, verbose = 0):
    try:
        assert not recompute
        from mmgroup.dev.mm_reduce.v1_mod3_data import V1_MOD3
    except (ImportError, AssertionError):
        result, _ = make_v1_mod3(verbose)
        write_v1_mod3(result)
        time.sleep(0.3)
        from mmgroup.dev.mm_reduce.v1_mod3_data import V1_MOD3
        
    if verbose:
        print("Cocode indices:", make_cocode_indices())
        v1 = MMV3('S',V1_MOD3) 
        print("A part of vector v1 (mod 3):")
        print(v1['A'])
        print("Weight of v1 is", len(V1_MOD3))
    assert check_v1_mod3(MMV3('S',V1_MOD3))
    return V1_MOD3

 

#######################################################################
# Compute tables for dealing with vevtor v1 (mod 3)
#######################################################################

def map_y(y_index):
    i, j = (y_index >> 14) & 0x1f, (y_index >> 8) & 0x1f
    vect = (1 << i) + (1 << j)
    gc = vect_to_cocode(vect)
    assert 0 <= gc < 0x800
    return gc 
    
    
def map_x(x_index):
    v2 = mm_aux_index_sparse_to_leech2(x_index) 
    return ((v2 & 0xfff) << 12) | ((v2 >> 12) & 0xfff)    


def make_v1_mod3_tags(V1_MOD3):
    ov = MMV3('S',V1_MOD3).data 
    tags_y = np.array(V1_MOD3[:11], dtype = np.uint32) 
    tags_x = np.array(V1_MOD3[11:11+24], dtype = np.uint32)
    watermark_perm = np.zeros(9, dtype = np.uint32)
    ok = mm_op_watermark_A(3, ov, watermark_perm)
    assert ok >= 0, ok

    solve_yt = np.zeros(11, dtype = np.uint64)
    assert len(tags_y) == 11
    nrows = 0
    for y in tags_y:
        eqn = map_y(y)
        nrows += leech2matrix_add_eqn(solve_yt, nrows, 11, eqn)
    assert nrows == 11, nrows
    solve_y = list(bitmatrix64_t(solve_yt, 11))
    assert len(solve_y) == 11
    assert mm_aux_mmv_extract_sparse_signs(3, ov, tags_y, 11) == 0
    
    solve_xt = np.zeros(24, dtype = np.uint64)
    assert len(tags_x) == 24
    nrows = 0
    for i, x in enumerate(tags_x):
        eqn =  map_x(x)  
        nrows += leech2matrix_add_eqn(solve_xt, nrows, 24, eqn)
    assert nrows == 24, nrows
    solve_x = list(bitmatrix64_t(solve_xt, 24))
    assert mm_aux_mmv_extract_sparse_signs(3, ov, tags_x, 24) == 0
    
    # Concatenate computed lists to the numpy array 'v1_mod3_tags'
    tags_x_leech2 = [mm_aux_index_sparse_to_leech2(x) for x in tags_x]
    assert 0 < min(tags_x_leech2)
    v1_mod3_tags = np.array(sum(map(list, [
        watermark_perm, tags_y, solve_y, tags_x_leech2, solve_x
    ]), []), dtype = np.uint32)
    assert len(v1_mod3_tags) == 9 + 2*(11+24), len(v1_mod3_tags)
    return v1_mod3_tags


#######################################################################
# Class containing tables for dealing with vevtor v1 (mod 3)
#######################################################################


class V1_Mod3_Table:
    directives = {}
    tables = None

    def compute_data(self):
        """Compute data for table classes

        We cannot to this ab initio; because the required inputs 
        are not available if sys.arg contains the argument ``mockup``.
        """
        if not self.__class__.tables is None:
            return
        if import_pending:
            import_all()
        v1_data =  get_v1_mod3_data()
        self.__class__.tables = {
            "V1_MOD3_DATA": v1_data,
            "V1_MOD3_TAG_DATA": make_v1_mod3_tags(v1_data)
        }

    def __init__(self, *args, **kwds):
        self.compute_data()


class Mockup_V1_Mod3_Table:
    a_ov =  np.array([0], dtype = np.uint32)
    directives = {}
    tables = {
        "V1_MOD3_DATA": a_ov,
        "V1_MOD3_TAG_DATA": a_ov
    }
    def __init__(self, *args, **kwds):
        pass


Tables = V1_Mod3_Table
MockupTables = Mockup_V1_Mod3_Table




#######################################################################
# Testing
#######################################################################

def make_test_samples(ntests):
    for i in range(ntests):
        yield MM0('r', 'G_x0')



OFS3_WATERMARK_PERM = 0
OFS3_TAGS_Y = 9
OFS3_SOLVE_Y = 20
OFS3_TAGS_X = 31
OFS3_SOLVE_X = 55


OFS3_Z = (116416 >> 5)


def do_check_sample(V1_MOD3, v):
    """Check if v  == MMV3('a', V1_MOD3)

    This destroys v.
    """
    v1_neg = np.array(V1_MOD3.copy(), dtype = np.uint32)
    for i in range(len(v1_neg)):
        v1_neg[i] ^= 3 
    mm_aux_mmv_add_sparse(3, v1_neg, len(v1_neg), v.data)
    assert mm_op_checkzero(3, v.data) == 0


def do_test_sample(v1, v1_mod3_tags, g, test_C, verbose = 0):
    w03 = mm_op_eval_A_rank_mod3(3, v1.data, 0)
    assert w03 >> 48 == 23 and w03 & 0xffffffffffff != 0, hex(w03)
    g1 = np.zeros(11, dtype = np.uint32)
    len_g1 = 0
    v = v1.copy() * g
    w3 = mm_op_eval_A_rank_mod3(3, v.data, 0)
    assert w3 >> 48 == 23 and w3 & 0xffffffffffff != 0, hex(w3)
    w_type4 = gen_leech3to2_type4(w3)
    assert w_type4 > 0, w_type4 
    len_g1 = gen_leech2_reduce_type4(w_type4, g1);
    if (len_g1 < 0):
        err = "gen_leech2_reduce_type4 failed with error %d"
        raise ValueError(err % len_g1)
    work_A = v.data[:24].copy()
    res = mm_op_word_tag_A(3, work_A, g1, len_g1, 1)
    if verbose:
        print("g1 =", MM0('a', g1[:len_g1]))
        print(str_watermark(work_A))
    assert res >= 0, res
    perm_num = mm_op_watermark_A_perm_num( 
        3, v1_mod3_tags[OFS3_WATERMARK_PERM:], work_A
    )
    assert perm_num >= 0, perm_num
    if (perm_num):
        g1[len_g1] = 0xA0000000 + perm_num;
        res = mm_op_word_tag_A(3, work_A, g1[len_g1:], 1, 1);
        assert res >= 0,  res;
        len_g1 += 1;
    if verbose:
        print("g product intermediate =", 
            Xsp2_Co1(g * MM0('a',g1[:len_g1])))

    v_y = mm_aux_mmv_extract_sparse_signs(3, work_A, 
        v1_mod3_tags[OFS3_TAGS_Y:], 11);
    assert v_y >= 0, v_y
    y = leech2matrix_solve_eqn(v1_mod3_tags[OFS3_SOLVE_Y:], 11, v_y);
    if (y > 0):
        g1[len_g1] = 0xC0000000 + y;
        res = mm_op_word_tag_A(3, work_A, g1[len_g1:], 1, 1);
        assert res >= 0, res;
        len_g1 += 1;
    g_x = Xsp2_Co1(g * MM0('a',g1[:len_g1]))
    if verbose:
        print("g product intermediate =", g_x)
    #g_x *= Xsp2_Co1('p', 1) # This generates an error in next statemnt
    g_x.as_Q_x0_atom()  # check that g_x is in Q_x0

    # Reduction in Q_x0 yet to be done!!!!!!!!

    gi = Xsp2_Co1('a',g1[:len_g1]) ** -1

    v_x = mm_aux_mmv_extract_x_signs(3, v.data, gi._data, 
        v1_mod3_tags[OFS3_TAGS_X:], 24) 
    assert v_x >= 0, v_x
    x = leech2matrix_solve_eqn(v1_mod3_tags[OFS3_SOLVE_X:], 24, v_x)
    v_sign = ((x >> 12) & 0x7ff) ^ (x & 0x800);
    aa = xsp2co1_elem_read_mod3(v.data[OFS3_Z:], gi._data, v_sign, 24) 
    sign = aa - 1
    assert 0 <= sign <= 1
    sign ^= uint64_parity(x & (x >> 12) & 0x7ff);
    x ^= sign << 24

    x ^= ploop_theta(x >> 12);
    if x & 0xfff:
        g1[len_g1] = 0x90000000 + (x & 0xfff); 
        len_g1 += 1
    xx = (x >> 12) & 0x1fff;
    if xx:
        g1[len_g1] = 0xB0000000 + xx;
        len_g1 += 1
    g_sign = Xsp2_Co1(g * MM0('a',g1[:len_g1]))
    if verbose:
        print("g product final =", g_sign)
    assert Xsp2_Co1(g * MM0('a',g1[:len_g1])) == Xsp2_Co1()
    g_inv = Xsp2_Co1('a', g1[:len_g1])
    if test_C:
       from mmgroup.mm_reduce import mm_order_find_Gx0_via_v1_mod3
       a_gi = np.zeros(10, dtype = np.uint32)
       l = mm_order_find_Gx0_via_v1_mod3(v.data, a_gi)
       assert l >= 0, hex(l)
       g_inv_C = Xsp2_Co1('a', a_gi[:l])
       assert g_inv_C == g_inv
    return g_inv


    
    
def test_samples(V1_MOD3, v1_mod3_tags, ntests=1, test_C=0, verbose=0):
    v1 =  MMV3('S',V1_MOD3)
    for i, g in enumerate(make_test_samples(ntests)):
        if verbose:
            print("Test %d" % (i+1))
            print("g =", g)
        gi = do_test_sample(v1, v1_mod3_tags, g, test_C, verbose)
        #do_check_sample(V1_MOD3, v1 * g * gi)


#######################################################################
# Main program
#######################################################################


def parse_args():
    from optparse import OptionParser
    parser = OptionParser()
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option("-r",  dest="recompute", action="store_true",
        help="Recompute vector v1 (mod 3)" )
    parser.add_option("-c",  dest="test_C", action="store_true",
        help="When testing: test also the C function" )
    parser.add_option("-t", type="int", dest="ntests",  metavar="NTESTS",
        default=0,
        help = "Test vector v1 with  NTESTS random samples in group G_x0" )
    parser.add_option("-v",  dest="verbose", action="store_true",
        help="Verbose operation" )
    
    options, args = parser.parse_args()
    return options, args




if __name__ == "__main__":
   import_all()
   options, args = parse_args()
   if options.recompute and options.test_C:
       err = """Options -r and -c are imcompatible. You must compile the
C file after recomputing with option -r!"""
       raise ValueError(err)
   V1_MOD3 = get_v1_mod3_data(options.recompute, options.verbose)
   v1_mod3_tags = make_v1_mod3_tags(V1_MOD3)
   table = V1_Mod3_Table()
   mockup_table = Mockup_V1_Mod3_Table()
   ntests = int(max(0, options.ntests))
   test_samples(V1_MOD3, v1_mod3_tags, ntests, options.test_C, 
       options.verbose)
   if ntests:
       print("%d tests passed" % ntests)


   



