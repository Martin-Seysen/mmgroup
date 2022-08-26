import os
import sys
import time
import numpy as np
import re
from collections import OrderedDict

from mmgroup.mat24 import vect_to_cocode
from mmgroup.mat24 import ploop_theta
from mmgroup.mm import mm_aux_index_sparse_to_leech2
from mmgroup.mm import mm_vector
from mmgroup.mm import mm_aux_mmv_extract_sparse_signs
from mmgroup.mm import mm_aux_mmv_extract_sparse
from mmgroup.structures.mm0_group import MM0
from mmgroup.mm_space import MMV, MMVector, MMSpace
from mmgroup.generators import mm_group_check_word_n
from mmgroup.generators import mm_group_words_equ
from mmgroup.generators import mm_group_n_mul_element
from mmgroup.generators import mm_group_n_reduce_element 
from mmgroup.generators import gen_leech3to2_type4
from mmgroup.generators import gen_leech2_reduce_type4
from mmgroup.clifford12 import chk_qstate12
from mmgroup.clifford12 import uint64_parity
from mmgroup.clifford12 import leech2matrix_add_eqn
from mmgroup.clifford12 import leech2matrix_solve_eqn
from mmgroup.clifford12 import bitmatrix64_t
from mmgroup.clifford12 import xsp2co1_half_order_word
from mmgroup.clifford12 import xsp2co1_power_word
from mmgroup.mm15 import op_eval_A_rank_mod3 as mm_op15_eval_A_rank_mod3 

from mmgroup.mm15 import op_copy as mm_op15_copy
from mmgroup.mm15 import op_compare as mm_op15_compare
from mmgroup.mm15 import op_word as mm_op15_word
from mmgroup.mm15 import op_norm_A as mm_op15_norm_A
from mmgroup.mm15 import op_watermark_A as mm_op15_watermark_A
from mmgroup.mm3 import op_watermark_A as mm_op3_watermark_A

from mmgroup.dev.mm_reduce import find_order_vector
from mmgroup.dev.mm_reduce.find_order_vector import stabilizer_vector
from mmgroup.dev.mm_reduce.find_order_vector import find_vector_71_mod3
from mmgroup.dev.mm_reduce.find_order_vector import find_vector_v94_mod5


MMV3 = MMV(3)
MMV15 = MMV(15)
MM = MM0  #  TODO: Fixme

ORDER_VECTOR = None
ORDER_TAGS = None


OFS_NORM_A = 0
OFS_DIAG_VA = 1
OFS_WATERMARK_PERM = 2
OFS_TAGS_Y = 26
OFS_SOLVE_Y = 37
OFS_TAGS_X = 48
OFS_SOLVE_X = 72
OFS_TAG_SIGN = 96


_DIR = os.path.split(find_order_vector.__file__)[0]
PY_FILENAME = os.path.join(_DIR, "order_vector_data.py")




#######################################################################
# Parse string to vector or group element
#######################################################################

m_vect = re.compile("MV<([0-9]+)")

def parse_str(value):
    if isinstance(value, str):
        mo = m_vect.match(value)
        if mo:
            p = int(mo.groups()[0])
            return MMVector(p, value)
        elif value.startswith("M0"):
            return MM0(value)
    return s



#######################################################################
# Assemble a test vector mod 15 from the input data
#######################################################################


def mm_element_from_data(data):
    if isinstance(data, str):
        result = parse_str(data)
        assert type(result) == MM0
        return result
    elif isinstance(data, MM0):
        return MM0(data)
    return MM0('a', data)

def vector_from_data(factor, data):
    if isinstance(data, str):
        v = parse_str(data)
    elif isinstance(data, MMVector):
        v = data
    else:
        v = MMV15('S', data)
    assert type(v) == MMVector, type(v)
    return MMV15(factor, v)

def identity_vector(factor):
    return factor * MMV15("+".join(["D_%s" % i for i in range(24)]))


   
def make_order_vector_mod3(s_g71, s_v71, s_gA, diag):
    v71 = vector_from_data(10, s_v71)
    g71 = mm_element_from_data(s_g71)
    w71 = stabilizer_vector(v71, g71, 71)
    assert w71 is not None
    w71 *= mm_element_from_data(s_gA)
    if not isinstance(diag, int): diag = diag[0]
    w71 += identity_vector(5 * diag % 15)
    diag = 0
    v3 = mm_op15_eval_A_rank_mod3(w71.data, diag) & 0xffffffffffff
    assert v3 != 0
    v_type4 = gen_leech3to2_type4(v3)
    assert v_type4 == 0x800000
    w71.reduce()
    return w71


def make_order_vector_mod15(w71, s_g94, s_v94):
    v94 = vector_from_data(6, s_v94)
    g94 = mm_element_from_data(s_g94)
    w94 = stabilizer_vector(v94 - v94 * g94, g94**2, 47)
    assert w94 is not None
    w = w71 + w94
    w.reduce()
    return w


def make_order_vector(s_g71, s_v71, s_gA, diag, s_g94, s_v94):
    v3 = make_order_vector_mod3(s_g71, s_v71, s_gA, diag)
    return make_order_vector_mod15(v3, s_g94, s_v94)





def map_y(y_index):
    i, j = (y_index >> 14) & 0x1f, (y_index >> 8) & 0x1f
    vect = (1 << i) + (1 << j)
    gc = vect_to_cocode(vect)
    assert 0 <= gc < 0x800
    return gc 
    
    
def map_x(x_index):
    v2 = mm_aux_index_sparse_to_leech2(x_index) 
    return ((v2 & 0xfff) << 12) | ((v2 >> 12) & 0xfff)    


def make_order_tags(order_vector, tags_y, tags_x, tag_sign):
    ov = order_vector.data
    tags_y = np.array(tags_y, dtype = np.uint32) 
    tags_x = np.array(tags_x, dtype = np.uint32)
    tag_sign =  np.array(tag_sign, dtype = np.uint32) 
    watermark_perm = np.zeros(24, dtype = np.uint32)
    ok = mm_op15_watermark_A(ov, watermark_perm)
    assert ok >= 0

    solve_yt = np.zeros(11, dtype = np.uint64)
    assert len(tags_y) == 11
    nrows = 0
    for y in tags_y:
        eqn = map_y(y)
        nrows += leech2matrix_add_eqn(solve_yt, nrows, 11, eqn)
    assert nrows == 11, nrows
    solve_y = list(bitmatrix64_t(solve_yt, 11))
    assert len(solve_y) == 11
    assert mm_aux_mmv_extract_sparse_signs(15, ov, tags_y, 11) == 0
    
    solve_xt = np.zeros(24, dtype = np.uint64)
    assert len(tags_x) == 24
    nrows = 0
    for i, x in enumerate(tags_x):
        eqn =  map_x(x)  
        nrows += leech2matrix_add_eqn(solve_xt, nrows, 24, eqn)
    assert nrows == 24, nrows
    solve_x = list(bitmatrix64_t(solve_xt, 24))
    
    # Concatenate computed lists to the global numpy array 'ORDER_TAGS'
    order_tags = np.array(sum(map(list, [
        [mm_op15_norm_A(ov), 0], 
        watermark_perm, tags_y, solve_y, tags_x, solve_x, tag_sign
    ]), []), dtype = np.uint32)
    assert len(order_tags) == 97, len(order_tags)
    t0 = mm_aux_mmv_extract_sparse_signs(
        15, ov, order_tags[OFS_TAGS_X:], 24)
    assert t0 == 0
    return order_tags



def order_vector_from_dict(d):
    S_G71, S_V71 = d["S_G71"], d["S_V71"]
    S_GA, DIAG_VA = d["S_GA"], d["DIAG_VA"]
    S_G94, S_V94 = d["S_G94"], d["S_V94"]
    data1 = [S_G71, S_V71, S_GA, DIAG_VA, S_G94, S_V94]
    ov = make_order_vector(*data1)
    TAGS_Y, TAGS_X =  d["TAGS_Y"], d["TAGS_X"]
    TAG_SIGN = d["TAG_SIGN"]
    data2 = [TAGS_Y, TAGS_X, TAG_SIGN]
    o_tag = make_order_tags(ov, *data2)
    keys =  [
      "S_G71", "S_V71", "S_GA", "DIAG_VA", "S_G94", "S_V94",
      "TAGS_Y", "TAGS_X", "TAG_SIGN"
    ]
    d1 = {k: d[k] for k in keys}
    return ov, o_tag, d1



#######################################################################
# Check that the test vector supports reduction
#######################################################################

Y_INDICES = [("A", i, j) for i in range(2) for j in range(i+1, 24)]
X_INDICES = [("B", i, j) for i in range(2) for j in range(i+1, 24)]
X_INDICES += [("C", 0, j)  for j in range(1, 24)]
BASIS = [0] + [1 << i for i in range(11)]
X_INDICES += [("X", i, j) for j in range(0, 24)  for i in  BASIS]
del BASIS


    
def eqn_system(vector, tag_indices, map_f, n):
    entries = [vector[index] for index in tag_indices]
    indices = [MMSpace.index_to_sparse(*x) for x in tag_indices]
    matrix= np.zeros(24, dtype = np.uint64)
    rows, cols = 0, n
    out_indices = []
    for (entry, index) in zip(entries, indices):
        if 1 <= entry < vector.p:
            eqn = map_f(index)
            new_rows = leech2matrix_add_eqn(matrix, rows, cols, eqn)
            if new_rows:
                out_indices.append(index)
                rows += new_rows
                if rows == n:
                    break
    if rows < n:
        return None
    mask = (1 << n) - 1
    out_matrix = [int(matrix[i]) & mask for i in range(n)] 
    return  out_indices   


def eqn_sign(vector):
    p = vector.p
    if p == 15:
        for j in range(24):
            if 1 <= vector["Z", 0, j] < 15:
                return [ MMSpace.index_to_sparse("Z", 0, j) ]
    if p == 3:
        for i in range(100):
            v2, v3 = vector["Z", i, 2:4]
            if 1 <= v2 < 3 and 1 <= (v2 - v3) % 3 < 3:
                return [ MMSpace.index_to_sparse("Z", i, 2) ]
    return None


def augment_v_data(v, data):
    a = np.array(data, dtype = np.uint32)
    mm_aux_mmv_extract_sparse(v.p, v.data, a, len(a))
    return [x for x in a]


def check_v(v, verbose = 0):
    vB = v['B']
    mark = np.zeros(24, dtype = np.uint32)
    if v.p == 15:
        for p in (3,5):
            if (vB % p == 0).all():
                if verbose:
                    print("Vector may be zero (mod %d)" % p)
                return None
        if mm_op15_watermark_A(v.data, mark) < 0:
            if verbose:
                print("Permutation watermarking failed")
            return None
    elif v.p == 3:
        if (vB == 0).all():
            if verbose:
               print("Vector may be zero (mod 3)")
            return None
        if mm_op3_watermark_A(v.data, mark) < 0:
            if verbose:
                print("Permutation watermarking failed")
            return None
    else:
        if verbose:
            print("Vector is not given mod 3 or mod 15")
        return None

    result_y = eqn_system(v, Y_INDICES, map_y, 11)
    if result_y is None:
        if verbose:
            print("Check Y failed")
        return result
    result_x = eqn_system(v, X_INDICES, map_x, 24)
    if result_x is None:
        if verbose:
            print("Check X failed")
        return result
    result_sign = eqn_sign(v)
    if result_sign is None:
        if verbose:
            print("Check for sign failed")
        return result
    results = [result_y, result_x, result_sign]
    return tuple([augment_v_data(v, data) for data in results])
    



#######################################################################
# Seach for the relevant data
#######################################################################

m_vect = re.compile("MV<([0-9]+)")

def to_list(value):
    if isinstance(value, str):
        value = parse_str(value)
    if isinstance(value, MM0):
        return list(map(int, value.mmdata))
    if isinstance(value, MMVector):
        return list(map(int, value.as_sparse()))
    if isinstance(value, int):
        return [value]
    if isinstance(value, list):
        return list(map(int, value))
    err = "Cannot convert type-'%s' object to list"
    raise TypeError(err % type(value))



def str_data(text, data):
    data = to_list(data)
    s = "%s = [\n   " % text
    for i, x in enumerate(data):
        s += hex(x) + ","
        s += "\n   " if i % 6 == 5 else " "
    s += "\n]\n"
    return s


def find_order_vector(verbose = 0):
    verbose = 1
    if verbose:
        print("Trying to find a vector of order 71")
    for trials in range(200,-1,-1):
        s_g71, s_v71, s_gA, diag = find_vector_71_mod3(verbose)
        v3 = make_order_vector_mod3(s_g71, s_v71, s_gA, diag)
        tag_data_mod3 = check_v(v3 % 3, verbose=verbose)
        if tag_data_mod3 is not None:
            if verbose:
                print("Tests for v71 passed")
            break
    if not trials:
        err = "No suitable vector mod 3 in the monster representation found"
        raise ValueError(err) 
    if verbose:
        print("Trying to find a vector of order 94")
    for trials in range(200,-1,-1):
        s_g94, s_v94 = find_vector_v94_mod5(verbose = verbose)
        if s_v94 is None:
            continue
        v = make_order_vector_mod15(v3, s_g94, s_v94)
        tag_data = check_v(v, verbose=verbose)
        if tag_data is not None:
            if verbose:
                print("Tests for v94 passed")
            break
    if not trials:
        err = "No suitable vector in the monster representation found"
        raise ValueError(err) 
    v_data =  s_g71, s_v71, s_gA, diag, s_g94, s_v94
    V_NAMES =  ["S_G71", "S_V71", "S_GA", "DIAG_VA", "S_G94", "S_V94"]
    TAG_NAMES =  [ "TAGS_Y", "TAGS_X", "TAG_SIGN"]   
    if verbose:        
        for text, data in zip(TAG_NAMES, tag_data):
            print(str_data(text, data))
    result = OrderedDict(zip(V_NAMES + TAG_NAMES, v_data + tag_data))
    return result


HEADER = """# This file has been created automatically, do not change!
# For documentation see module mmgroup.dev.mm_reduce.order_vector.

"""

def write_order_vector(result):
    print("Writing file " + PY_FILENAME)
    f = open(PY_FILENAME, "wt")
    print(HEADER, file = f)
    for text, data in result.items():
        print(str_data(text, data), file = f)
    f.close()
    



def get_order_vector(recompute = False, verbose = 0):
    try:
        from mmgroup.dev.mm_reduce import order_vector_data
        assert not recompute
    except (ImportError, AssertionError):
        result = find_order_vector(verbose)
        write_order_vector(result)
        from mmgroup.dev.mm_reduce import order_vector_data
    ov, o_tag, data = order_vector_from_dict(order_vector_data.__dict__)
    return ov, o_tag, data






