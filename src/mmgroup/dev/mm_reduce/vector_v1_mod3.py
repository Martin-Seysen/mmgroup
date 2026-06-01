from collections import defaultdict
import numpy as np

if __name__ == "__main__":
    import sys
    import os
    sys.path.append(os.path.join('..', '..', '..'))    

try:
    import mmgroup
    from mmgroup import MMV, MM0, MM, MMSpace, Xsp2_Co1, mat24
    from mmgroup.generators import gen_leech2_reduce_type4
    from mmgroup.generators import gen_leech2_op_word_many
    from mmgroup.generators import gen_leech3to2_type4
    from mmgroup.clifford12 import leech2matrix_add_eqn
    from mmgroup.clifford12 import bitmatrix64_t
    from mmgroup.clifford12 import leech2matrix_solve_eqn
    from mmgroup.clifford12 import xsp2co1_rep_mod3_unit_vector
    from mmgroup.clifford12 import xsp2co1_rep_mod3_mul_word
    from mmgroup.clifford12 import xsp2co1_rep_mod3_find_nonzero
    from mmgroup.mm_op import mm_aux_index_sparse_to_leech2
    from mmgroup.mm_op import mm_op_eval_A_rank_mod3
    from mmgroup.mm_op import mm_aux_mmv_add_sparse
    from mmgroup.mm_op import mm_op_word_tag_A
    from mmgroup.mm_op import mm_aux_index_leech2_to_sparse
    from mmgroup.mm_op import mm_aux_mmv_extract_sparse
    mmgroup_present = True
except (ImportError, ModuleNotFoundError):
    mmgroup_present = False

mmgroup_reduce_present = False
if mmgroup_present:
    try:
        from mmgroup.mm_reduce import mm_order_find_Gx0_via_v1_mod3
        mmgroup_reduce_present = True
    except:
        pass



####################################################################
####################################################################
# Creating the vector v_1 (mod 3)
####################################################################
####################################################################






####################################################################
# Part A of vector v_1
####################################################################


def make_A12():
    """Return a certain symmetric 12 times 12 matrix A12

    Matrix A12 satisfies the following properties:

    All entries of A12 are either 0 or 1.

    For row i of matrix A let h(i, A) be the pair (A[i,i], w(i)),
    where w(i) is the weight of row i, i.e. the number of its
    nonzero entries. Then we have h(i, A12) != h(j, A12) for i != j.

    The determinant of A12 is visibly equal to one. To check this,
    decompose A12 into four blocks of 6 times 6 matrices.
    """
    a = np.zeros((12, 12), dtype = np.int32)
    for i in range(12):
        a[i,i] = 1 if i < 6 else 0
        for j in range(i + 6, 12):
            a[j, i] = a[i, j] = 1
    det = np.round(np.linalg.det(a))
    assert abs(1.0 - det) < 1.0e-8
    return a

A12 = make_A12()

    
MAP_A12 = [0,1,2,3,4,5,6,8,9,10,12,16]
KER_A12 = [7]

def data_A24():
    for i in range(12):
        ii = MAP_A12[i]
        for j in range(i, 12):
            jj, e = MAP_A12[j], A12[i,j]
            if e:
                 yield (int(e) % 3, 'A', ii, jj)
    neg_diag = set(MAP_A12) | set(KER_A12)
    for i in range(24):
        if i not in neg_diag:
            yield (2, 'A', i, i)



def make_A24():
    A24 = np.zeros((24,24), dtype = np.uint32)
    for x, tag, i, j in data_A24():
        if tag == "A":
           A24[i,j] = A24[j,i] = x
    return A24


def process_A24(A24, verbose = 0):
    """Compute a certain 24 times 24 matrix

    The function returns a quadruple ``(A, sp, h, eval_v_y)``

    Here ``A`` is a symmetric 24 times 24 Matrix with entries 0, 1, 2
    to be interpreted as integers modulo 3. The kernel of ``A``
    (modulo 3) has dimension one and is spanned by a unit vector.
    The return value ``sp`` is the list of nonzero entries of ``A``
    in *sparse* notation.

    Matrix ``A`` is constructed from matrix ``A12`` given by function
    ``make_A12``. For row ``i`` of matrix ``A`` we define the pair
    ``h(i, A)`` as in that function. Then the return value ``h`` is a
    dictionary mapping the pair ``h(i, A)`` to ``i``, provided that
    ``i`` is  uniquely defined by ``h(i, A)``. This works for
    ``0 <= i <= 10`` and ``i = 12, 16``.

    The return value ``eval_v_y`` is a list of pairs ``(i, j)`` with
    ``0 <= j < i < 24`` such that  ``A[i, j] = 1``, and the Golay
    cocode elements corresponding the the pairs [i,j] in that
    list are a basis of the even part of the Golay cocode.
    The entries of the list are encoded as sparse vectors.
    """
    sp, eval_v_y = [], []
    for x, tag, i, j in data_A24():
        if tag == "A":
           sp.append(MMSpace.index_to_sparse(tag, i, j) + (x & 3))
           if i < j and (i == 0 or j == 16):
               assert x == 1
               eval_v_y.append(MMSpace.index_to_sparse('A', j, i))
    assert (A24[7] == 0).all()
    A23 = np.copy(A24)
    A23 = np.delete(A23, 7, 0)
    A23 = np.delete(A23, 7, 1)
    det = np.linalg.det(A23)
    idet = int(np.round(det))
    assert abs(det - idet) < 1.0e-5
    assert idet % 3 != 0
    h, n = {}, 0
    for i in range(24):
        if A24[i,i] < 2:
            hash = int(A24[i,i]), int(np.count_nonzero(A24[i]))
            assert hash not in h
            h[hash] = i
        else:
            n += 1
    assert set(h.values()) == set(MAP_A12) | set(KER_A12)
    assert len(h) + n == 24
    sp.sort()
    eval_v_y.sort()
    return sp, eval_v_y





def make_hash_a(A24):
    a = np.array(A24, dtype = np.int32) % 3
    n_hashes, ha = defaultdict(int), []
    for i in range(24):
        hash = 32 * int(a[i,i]) + int(np.count_nonzero(a[i]))
        ha.append(hash)
        n_hashes[hash] += 1
    HASHED = [0,1,2,3,4,5,8]
    hash_table = []
    for i, entry in enumerate(HASHED):
        hv = ha[entry]
        assert n_hashes[hv] == 1
        hash_table.append(hv)
    return ha, hash_table





####################################################################
####################################################################
# Compute data for obtaining the 2-subgroup 2^{1+24+11} of N_x0
####################################################################
####################################################################




####################################################################
# Compute y_d from part A of the vector v_1
####################################################################



def map_y(eval_v_y):
    solve_yt = np.zeros(11, dtype = np.uint64)
    assert len(eval_v_y) == 11
    nrows = 0
    eval_y = []
    for y in eval_v_y:
        i, j = (y >> 14) & 0x1f, (y >> 8) & 0x1f
        eqn = mat24.vect_to_cocode((1 << i) + (1 << j))
        s = leech2matrix_add_eqn(solve_yt, nrows, 11, eqn)
        assert s == 1
        nrows += 1
        eval_y.append([i,j])
    solve_y = list(int(x) for x in bitmatrix64_t(solve_yt, 11))
    assert len(solve_y) == 11
    eval_y = np.array(eval_y, dtype = np.uint8)
    return eval_y, np.array(solve_y, dtype = np.uint32)

####################################################################
#  Compute x_r from parts BCTX of the vector v_1
####################################################################


def leech2_covector(v2):
    return ((v2 & 0xfff) << 12) | ((v2 >> 12) & 0xfff)



def tuple_data_BCTX():
    for i in range(1,12):
        yield (1, 'B', MAP_A12[0], MAP_A12[i])
    yield (1, 'C', 2, 3)
    for i in range(11):
        o = mat24.gcode_to_octad(1 << i, 0)
        yield (1, 'T', o, 0)
    yield (1, 'X', 0, 0)


def data_BCTX():
    sp_x = []
    for x, tag, i, j in tuple_data_BCTX():
        sp = MMSpace.index_to_sparse(tag, i, j)
        assert sp > 0
        sp_x.append(sp + (x & 3))
    return sp_x


def map_x(eval_v_x):
    eval_x = []
    solve_xt = np.zeros(24, dtype = np.uint64)
    nrows = 0
    for v_x in eval_v_x:
        leech2 = mm_aux_index_sparse_to_leech2(v_x)
        eval_x.append(leech2)
        co_v2 = leech2_covector(leech2)
        assert leech2matrix_add_eqn(solve_xt, nrows, 24, co_v2) == 1
        nrows += 1
    assert len(eval_x) == 24
    assert nrows == 24, nrows
    eval_x = np.array(eval_x, dtype = np.uint32)
    solve_x = list(bitmatrix64_t(solve_xt, 24))
    return eval_x, np.array(solve_x, dtype = np.uint32)





####################################################################
# Compute sign from parts ZY of the vector v_1
####################################################################

def process_Z():
    x, tag, i, j = 1, 'Z', 0, 0
    sparse = MMSpace.index_to_sparse(tag, i, j)
    return [sparse + (x & 3)]




####################################################################
# obtain dict of constantans for compute 2-subgroup of Co_1
####################################################################

def get_eval_dict(sp):
    """Return a dictionary ``d`` for obtaining 2-subgroup of G_x0

    Input ``sp`` is a list of 11+24+1 coordinate positions of an
    order vector ``v1`` in the rep of the Monster mod 3. The order
    vector ``v1``  must have value 1 at all these positions.

    Let math:`Q_{x0}` be the subgroup  of structure math:`2^{1+24}`
    of math:`N_{x0}` (of structure math:`2^{1+24+11}.M_{24}`) as
    usual; and let math:`N_2` be the subgroup of math:`N_{x0}` of
    structure math:`2^{1+24+11}`. Let ``v1`` be an order vector and
    ``v`` the the order vector transformed by an element of math:`N_2`.

    The first 11 positions in ``sp`` must be at part 'A' of the
    vector. They must be independent in the sense that they uniquely
    determine an ``y_d``, with ``d`` in the Parker loop (modulo its
    centre),  such that ``y_d`` transforms ``v`` to a vector
    in ``v1`` * math:`Q_{x0}`. Output d["EVAL_Y"] is a list of
    11 pairs (i,j) corresponding to the positions of the entries in
    part 'A' given by the 11 postions in ``sp``. Output d["EQU_Y"] 
    is an 11 times 11 bit matrix. Multiplying the vector of the signs
    of the 11 entries described above by that matrix yields ``d``.

    in the sequel we assume that ``v`` is in ``v1`` * math:`Q_{x0}`.  
    The next 24 positions in ``sp`` must be at parts 'BCTX' of the
    vector. They must be independent in the sense that they uniquely
    determine an element ``r`` of the Leech lattice mod 2, such that
    ``x_r`` transforms ``v`` to a vector  ``v1`` * :math:`x_{\pm 1}`.
    Output d["EVAL_X"] contains these 24 position in Leech lattice
    encoding.  Output d["EQU_Y"]  is a 24 times 24 bit matrix.
    Multiplying the vector of the signs of the 24 entries described
    above by that matrix yields ``r``.

    Finally, d["SP_Z"][0] equal to the last entry in the list ``sp``.
    This is a coordinate position in party 'ZY' that can be used to
    obtain the sign of a ``v`` in  ``v1`` * :math:'x_{\pm 1}``.
    """
    d = {}
    d["EVAL_Y"], d["EQU_Y"]  =  map_y(sp[:11])
    d["EVAL_X"], d["EQU_X"]  =  map_x(sp[11:11+24])
    d["SP_Z"] = sp[35:36] 
    return d

def table_from_eval_dict(d):
    NAMES =  ["EVAL_Y", "EQU_Y", "EVAL_X", "EQU_X", "SP_Z"]
    return [d[x] for x in NAMES]
    

####################################################################
####################################################################
# Class ReduceGx0Data, for assembling the parts of the vector v_1
####################################################################
####################################################################





def display_part_A(cls):
    print("Data related to part 'A' of the Griess algebra")
    print("Submatrix A12 of part A of vector v_1:")
    print(cls.A12)
    print("Determinant of A12 is %.3f" % np.linalg.det(A12))

    print("Part A of vector v_1 (modulo 3):")
    print(cls.A24)
    A23 = np.copy(cls.A24)
    A23 = np.delete(A23, 7, 0)
    A23 = np.delete(A23, 7, 1)
    print("Determinant of image of that part is %.3f"
             % np.linalg.det(A23))
    if not mmgroup_present:
        s = "Cannot display more data, since mmgroup is not avialable"
        print(s)
        return

    print("Hash values for part A:")
    for i, (key, value) in enumerate(cls.HASH.items()):
        print("%s: %2d," % (key, value),
            end = "\n" if i % 5 == 4 else " ")
    print()

    print("Matrix SP in sparse notation (hex):")
    for i, e in enumerate(cls.SP):
        print("%07x" % e, end = "\n" if i % 8 == 7 else " ")
    print()


    print("Values of part A for determining generator y:")
    for i, d in enumerate(cls.EVAL_Y):
        print("%7s, " % (d,), end = "\n" if i % 6 == 5 else " ")
    print()
    print("Matrix for solving y (hex):")
    for d in cls.EQU_Y:
        print("%03x, " % d, end = "")
    print()
    print("Leech2 values of parts BCTX for finding generator x (hex):")
    for i, d in enumerate(cls.EVAL_X):
        print("%06x" % d, end = "\n" if i % 8 == 7 else " ")
    print("Matrix for solving x (hex):")
    for i, d in enumerate(cls.EQU_X):
        print("%06x" % d, end = "\n" if i % 8 == 7 else " ")
    print()





class ReduceGx0Data:
    A12 = A12
    A24 = make_A24()
    if mmgroup_present:
        SP_Y,  EVAL_V_Y = process_A24(A24)
        #EVAL_Y, EQU_Y = map_y(EVAL_V_Y)
        SP_X = data_BCTX()
        #EVAL_X, EQU_X = map_x(SP_X)
        SP_Z = process_Z()
        SP = SP_Y + SP_X + SP_Z
        SP = np.array([len(SP)] + SP, dtype = np.uint32)
        HASH, HASH_A = make_hash_a(A24)
        D =  get_eval_dict(EVAL_V_Y + SP_X + SP_Z)
        TABLE1 = [HASH, HASH_A]
        TABLE2 = table_from_eval_dict(D)

    @classmethod
    def display(cls):
        return display_part_A(cls)

    @classmethod
    def set_class_attr(cls, name, value):
        setattr(cls, name, value)



if 0 and mmgroup_present:
     d = ReduceGx0Data.D
     for name in ["EVAL_Y", "EQU_Y", "EVAL_X", "EQU_X", "SP_Z"]:
         ReduceGx0Data.set_class_attr(name, d[name])


class MockupReduceGx0Data:
    TABLE1 = [[0]]
    TABLE2 = [[0,0]]



####################################################################
####################################################################
# Tables for code generation
####################################################################
####################################################################

class Tables:
    directives = {}
    tables = None

    def __init__(self, *args, **kwds):
        """Compute data for table classes

        We cannot to this ab initio; because the required inputs 
        are not available if sys.arg contains the argument ``mockup``.
        """
        if not self.__class__.tables is None:
            return
        if not mmgroup_present:
            ERR = "Cannot create tables, since mmgroup is not present"
            raise ImportError(ERR)
        self.__class__.tables = {
            "V1_MOD3_VECTOR": ReduceGx0Data.SP,
            "V1_MOD3": ReduceGx0Data
        }



class MockupTables:
    a_ov =  np.array([0], dtype = np.uint32)
    directives = {}
    tables = {
        "V1_MOD3_VECTOR": [0],
        "V1_MOD3": MockupReduceGx0Data
    }
    def __init__(self, *args, **kwds):
        pass




####################################################################
####################################################################
# Using the vector v_1 (mod 3)
####################################################################
####################################################################

####################################################################
# State-based recuction function
####################################################################



def process_v_1(hash_a = None, d = None):
    if hash_a is None:
       hash_a = ReduceGx0Data.HASH_A
    if d is None:
       d = ReduceGx0Data.D
    state = 1
    len_data = 1
    data = np.zeros(24, dtype = np.uint32)
    g = np.zeros(12, dtype = np.uint32)
    len_g = 0
    yield state, len_data, data, g[:len_g]

    # Here data[0] should contain the kernel of part 'A' of vector
    # v_1 (which is of type 4) as a vector in the Leech lattice
    # mod 2 in *Leech lattice encoding*.
    g = np.zeros(12, dtype = np.uint32)
    len_g = gen_leech2_reduce_type4(data[0], g);
    if len_g < 0 or len_g > 6:
        ERR = "Kernel of part A of of vector v_1 is not of type 4"
        raise ValueError(ERR)
    np.copyto(data[:len_g], g[:len_g])
    state = 2
    len_data = len_g
    # Yield an element h of :math:`G_{x0}` that transforms the input
    # axis to an image of the standard axis under an element of
    # :math:`n_{x0}`. Element h is given as a product of generators
    # in the array data[:len_data].
    yield state, len_data, data, g[:len_g]

    # Let 'Ah' be part 'A' of the input axis transformed by h.
    # Here data[i], 0 <= i < 24 should contain the following
    # information about the rows of matrix 'Ah':
    # data[i] = 32 * Ah[i,i] + w(i, Ah).
    # Here w(i, Ah) is the weight of row i of Ah, i.e. the number
    # of its nonzero entries.  Ah[i,i] must be 0, 1, or 2.
    pi = np.zeros(24, dtype = np.uint8)
    pi.fill(24)
    for i in range(24):
        for j in range(7):
            if data[i] == hash_a[j]:
                pi[j] = i
    pi[8], pi[6] = pi[6], 24
    res = mat24.perm_complete_heptad(pi)
    assert res == 0
    mat24.perm_check(pi)
    state = 3
    len_data = len(d["EVAL_Y"])
    for n, (i, j) in enumerate(d["EVAL_Y"]):
        data[n] = 0x2000000 + (int(pi[i]) << 14) + (int(pi[j]) << 8)
    pi_num = mat24.perm_to_m24num(pi)
    if pi_num:
        g[len_g] = 0xA0000000 + pi_num
        len_g += 1
    # Yield a list of ``len_data`` entries to be read from
    # matrix 'Ah', with 'Ah' as above. Indices of these entries
    # are encoded in *sparse notation*.
    yield state, len_data, data, g[:len_g]

    # Here entry data[k] should have been replaced by the value
    # of the entry Ah (mod 3) at the index given by the input
    # data[k] returned at the previous yield statement. Input
    # data[k] is given in *sparse notation*, with k running
    # from 0 to data_len - 1.
    err, v_y = 2, 0
    for k in range(len_data):
        err &= data[k] + 1
        v_y |= (data[k] & 2) << k
    assert err == 2, err
    v_y >>= 1
    y = leech2matrix_solve_eqn(d["EQU_Y"], 11, v_y)
    assert y >= 0
    if y > 0:
        g[len_g] = 0xC0000000 + y
        len_g += 1

    # Now we are in the group Q_x0.
    for i in range(24):
        data[i] = d["EVAL_X"][i]
    g_inv = np.zeros(12, dtype = np.uint32)
    for i in range(len_g):
        g_inv[i] = g[len_g - 1 - i] ^ 0x80000000
    res = gen_leech2_op_word_many(data, 24, g_inv, len_g)
    assert res == len_g
    signs = sum(((x >> 24) & 1) << i for i, x in enumerate(data[:24]))
    for i in range(24):
        data[i]=  mm_aux_index_leech2_to_sparse(data[i]) & 0xffffff00
        assert data[i] > 0
    state = 4
    len_data = 24
    # Yield a list of ``len_data`` entries of the input axis
    # in *sparse notation* .
    yield state, len_data, data, g[:len_g]

    # Here entry data[k] should have been replaced by the value
    # of the entry of the input axis (mod 3) at the index given
    # by the input data[k]. These input data are given
    # in *sparse notation*, with k running
    err, v_x = 2, 0
    for k in range(24):
        err &= data[k] + 1
        v_x |= (data[k] & 2) << k
    assert err == 2, err
    v_x >>= 1
    v_x ^= signs
    x = leech2matrix_solve_eqn(d["EQU_X"], 24, v_x)
    x ^= mat24.ploop_theta((x >> 12) & 0xfff)
    if (x & 0xfff):
        g[len_g] = 0x90000000 + (x & 0xfff)
        len_g += 1
    x = (x >> 12) & 0x1fff
    g[len_g] = 0xB0000000 + x
    len_g += 1
    rep = np.zeros(14, dtype = np.uint64)
    v_z = d["SP_Z"][0]
    i1, j1 = (v_z >> 14) & 0xfff, (v_z >> 8) & 0x1f
    sign = (v_z & 2) >> 1
    xsp2co1_rep_mod3_unit_vector(i1, 1 << (j1 + 32*sign), rep)
    for i in range(len_g):
        g_inv[i] = g[len_g - 1 - i] ^ 0x80000000
    res = xsp2co1_rep_mod3_mul_word(rep, g_inv, len_g)
    assert res >= 0
    # Now everythng is fine up to sign
    res = xsp2co1_rep_mod3_find_nonzero(rep);
    assert res >= 0
    signs = (res >> 1) & 1
    data[0] = res & 0xffffff00
    state = 5
    len_data = 1
    yield state, len_data, data, g[:len_g]

    # Same action expected as at pevious state 4
    assert (data[0] + 1) & 2 == 2
    signs = (signs ^ (data[0] >> 1)) & 1
    g[len_g-1] ^= signs << 12
    if (g[len_g-1] & 0x1fff) == 0:
        len_g -= 1
    # Finally, we will return the result g
    for i in range(len_g):
       data[i] = g[i]
    state = 6
    len_data = len_g
    yield state, len_data, data, g[:len_g]
    # Return the final result g


####################################################################
# Using the State-based recuction function
####################################################################


def mul_g_data(g, data, verbose = 0):
    h = g.__class__('a', data)
    gh = g.copy() * h
    if verbose:
        print([hex(x) for x in gh.mmdata])
    return gh

def watermark_row(a, i):
    a = int(a)  # numpy problem
    d = (a >> (2*i)) & 3
    s = (a ^ (a >> 1)) & 0x555555555555
    s = (s + (s >> 2)) & 0x333333333333
    s = (s + (s >> 4)) & 0x0F0F0F0F0F0F
    s += s >> 8
    s += s >> 16
    s += s >> 32
    return (d << 5) + (s & 0x1f)

def py_reduce_v_1(v3, g = None, verbose = 0):
    st = process_v_1()
    state, length, data, h = next(st)
    w3 = mm_op_eval_A_rank_mod3(3, v3.data, 0)
    data[0] = gen_leech3to2_type4(w3)
    assert data[0] > 0
    state, length, data, h = next(st)
    assert mul_g_data(g, h, verbose).subtype == (4, 8)
    A = np.array(v3.data[:24], dtype = np.uint64)
    res = mm_op_word_tag_A(3, A, data, length, 1)
    assert res >=  0
    for i, a in enumerate(A):
        data[i] =  watermark_row(A[i], i)
    state, length, data, h = next(st)
    mul_g_data(g, h, verbose)
    mm_aux_mmv_extract_sparse(3, A, data, length)
    state, length, data, h = next(st)
    assert state == 4
    mul_g_data(g, h, verbose)
    mm_aux_mmv_extract_sparse(3, v3.data, data, length)
    state, length, data, h = next(st)
    assert state == 5
    mul_g_data(g, h, verbose)
    mm_aux_mmv_extract_sparse(3, v3.data, data, length)
    state, length, data, h = next(st)
    assert state == 6
    id = mul_g_data(g, h, verbose)
    assert id == id.__class__()
    return h

def C_reduce_v_1(v3):
    h = np.zeros(12, dtype = np.uint32)
    res = mm_order_find_Gx0_via_v1_mod3(v3.data,  h)
    assert res >= 0, (res, hex(res))
    return h[:res]

####################################################################
####################################################################
# Creating the vector v_1 (mod 3)
####################################################################
####################################################################



def py_v_1_mod3():
    v3 = MMV(3)()
    mm_aux_mmv_add_sparse(3, ReduceGx0Data.SP[1:],
        ReduceGx0Data.SP[0], v3.data)
    return v3


####################################################################
####################################################################
# Test reduction in G_x0
####################################################################
####################################################################

####################################################################
# Test python version
####################################################################


W_no_C = "C function mm_order_find_Gx0_via_v1_mod3 not implemented"


def one_py_test_reduce_Gx0(g, test_C = True, verbose = 0):
    global W_no_C
    v3 = py_v_1_mod3()
    v3 *= g
    h = py_reduce_v_1(v3, g)
    if test_C:
        if mmgroup_reduce_present:
            hc = C_reduce_v_1(v3)
            if not (h == hc).all():
                print("Expected:", [hex(x) for x in h])
                print("Obtained:", [hex(x) for x in hc])
                E = "Error in function mm_order_find_Gx0_via_v1_mod3"
                raise ValueError(E)
            elif verbose > 1:
                print("Test of C function passed")
        else:
            print(W_no_C)
    else:
        print("Test of C function skipped")

def py_test_Gx0(ntests = 50, test_C = True, verbose = 0):
    print("\nTest function mm_order_find_Gx0_via_v1_mod3()")
    if not mmgroup_present:
        s = "Cannot test reduction in G_x0, since mmgroup is not present"
        print(s)
        return
    for i in range(ntests):
        g = Xsp2_Co1('r', 'G_x0')
        if verbose:
            print("Test %d, g = " % (i+1), g)
        one_py_test_reduce_Gx0(g, test_C, verbose)



####################################################################
####################################################################
# Main program for testing
####################################################################
####################################################################



if __name__ == "__main__":
    ReduceGx0Data.display()
    py_test_Gx0(ntests = 50, test_C = 1, verbose = 2)


