from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import numpy as np
from numbers import Integral
from random import randint


import pytest


from mmgroup.structures.abstract_mm_rep_space import AbstractMmRepSpace
from mmgroup.structures.abstract_mm_rep_space import sparse_subspace
from mmgroup.mm_space import characteristics
from mmgroup.mm_op import mm_aux_index_sparse_to_extern
from mmgroup.mm_op import mm_aux_index_extern_to_sparse
from mmgroup.mm_op import mm_aux_mmv_extract_sparse
from mmgroup.mm_op import mm_op_word_tag_A

from mmgroup.tests.spaces.spaces import MMTestSpace
from mmgroup.tests.spaces.sparse_mm_space import SparseMmV
from mmgroup.tests.groups.mgroup_n import MGroupNWord


#from mm_vector import FastMmRepSpace, PRIMES
#from mm_vector import SparseMmRepVector
#from mm_vector import AbstractRepVector
#from mm_vector import mm_sparse_purge
#from mm_vector import basis_vectors_to_sparse

PRIMES = characteristics()


np.set_printoptions(formatter={'int':hex})

########################################################################
# group_blocks
########################################################################

"""Finding blocks in the matrix of an element g of the monster

Let R be the representation of the monster Mm and let g in Mm. We want 
to test the fast C implementation of the operation of g against the 
slow python implementation. While the C implementation always operates 
on a vector v in R of full length, the python implementation works on 
the individual components ov v, with run time depending on the number 
of componennts of v. So we select a uniform random vector v and an 
atom of the group g and we let the C implementation of g operate of v.
Then we let the python implementation of g operate on a few components
of v selected at random, and we compare the results. Of course, one
component of g(v) may depend on several components of v, and we
must apply the operation of g to (at least) all these components
for obtaning a correct component of g(v).

Let B be a (usually small) set of basis vectors of R and let g be 
an element of Mm. Then function group_blocks() maps B to a superset 
B_g of B of basis vectors of R with the following property:

g maps the linear span of B_g to a subspace spanned by basis vectors.

For a random vector v in R, let v1 be the the projection of v onto
the linear span of B_g. Then g(v1) is the projection of g(v) onto the 
space spanned by the basis vectors for which the coordinate in 
g(v1) is not zero.

If B_g is not too large, then g(v1) can also be computed with the 
python implementation without too much effort. This makes it
possible to test C implementation against the python implementation.

In the worst case, B_g may be the whole basis of R. Function 
group_blocks() raises ValueError if the computation of B_g is too
complicated. The function returns with success if the following 
conditions are met:

  -  g contains at most one non-monomial atom. Non-monomial 
     atoms are those with tag 't' or 'l'.
  -  An atom with tag 'l' has no direct or indirect predecessor 
     with tag 'p'.

"""


def group_blocks_t(a):
    """Auxiliary function for function group_blocks()

    Assuming that the numpy array 'a' aencodes a basis B, the 
    function returns the expanded basis B_g in case that an  
    atom with tag 't' has been found in the group element.

    The function is implemented as a generator. It yields
    unsigned 32-bit integers representing basis vectors.
    """
    ABC_tags = set()
    T_tags = set()
    for sp in a:
        tag = (sp >> 25) 
        if 1 <= tag <= 3:
            i, j = (sp >> 14) & 0x7ff, (sp >> 8) & 0x3f
            if i != j:
                ABC_tags.add(sp & 0x1ffff00)
            else:
                yield sp
        elif tag == 4:
            T_tags.add(sp & 0xffffc000 )
        else:
            yield sp
    for sp in ABC_tags:
        for tag in (1,2,3):
            yield (tag << 25) + sp
    for sp in T_tags:
        for j in range(64):
            yield sp + (j << 8)
            

def group_blocks_l(a):
    """Auxiliary function for function group_blocks()

    Assuming that the numpy array 'a' ancodes a basis B, the 
    function returns the expanded basis B_g in case that an  
    atom with tag 'l' has been found in the group element.

    The function is implemented as a generator. It yields
    unsigned 32-bit integers representing basis vectors.
    """
    A_tags = set()
    ZY_tags = set()
    for sp in a:
        tag = (sp >> 25) 
        if tag == 1:
            A_tags.add(sp & 0xffff1c00)
        elif tag in (6,7):
            ZY_tags.add(sp & 0xfffc1c00)
        else:
            yield sp
    for sp in A_tags:
        i, j = (sp >> 14) & 0x7ff, (sp >> 8) & 0x3f
        if i != j:
            for i0 in range(4):
                for j0 in range(4):
                    yield sp + (i0 << 14) + (j0 << 8)
        else:
            for i0 in range(4):
                for j0 in range(i0 + 1):
                    yield sp + (i0 << 14) + (j0 << 8)
    for sp in ZY_tags:
        for i0 in range(16):
            for j0 in range(4):
                yield sp + (i0 << 14) + (j0 << 8)




def group_blocks(basis_vectors, g):
    """Find block in matrix g operating on the basis_vectors.

    Here g in a group element of the monster group Mm, where the
    monster group is an instance of class MgroupN.

    'basis_vectors' represents a set B of basis vectors of the
    representatiom of Mm. It may be

     - a vector in a representation R_p of Mm, where R_p is an
       instance of (a subclass of) class AbstractMmRepSpace. 
       Then all basis vectors corresponding to nonzero components 
       of that vector are taken
     - A string of tags. Then for each tag a random basis
       vector with that tag is generated.
     - An array like object that will be converted to a one-
       dimensional numpy array of dtype np.uint32. This array
       is interpreted as a sparse representation of a vector v.
       In this case all basis vectors contained in the sparse 
       representation ar taken, ignoring coordinates.


    Given the set B of basis vectors and an element g of Mm, the 
    function returns a list B_g of basis vectors in sparse format
    as a numpy array of dtype = np.uin32.

    Here B_g is a superset the set B of basis vectors with the 
    following property:

    g maps the linear span of B_g to a subspace of the 
    representation of Mm spanned by basis vectors.

    The function may raise ValueError if the computation of B_g 
    is too complicated, see comment in the header of this 
    section for details. 
    """
    if isinstance(basis_vectors, str):
        a = sparse_subspace(*((x,) for x in basis_vectors))
    else:
        a = sparse_subspace(basis_vectors)
    # Now a is the list of basis vectors encoded as a numpy array.
            
    dirty_t = dirty_l = error = False
    # Here dirty_t (or dirty_l) means that we can no longer control
    # the proliferation of components due to coming atoms with 
    # tag 't' (or tag 'l').
    expansion_function = None

    for tag, i0 in g.as_tuples():
        if tag == 'p' and i0:
            dirty_l = True 
        elif  tag == 't' and i0 % 3:
            dirty_l = True 
            error = error or dirty_t
            expansion_function = group_blocks_t
        elif  tag == 'l' and i0 % 3: 
            dirty_t = True 
            error = error or dirty_l
            expansion_function = group_blocks_l
    if error:
        raise ValueError("Matrix too complicated for finding blocks")
    if expansion_function:
        a = np.fromiter(expansion_function(a), dtype = np.uint32)
    return a
        


########################################################################
# tests with a sparse vector v
########################################################################


def one_test_op(v, g, f_mul = None, verbose = 0):
    """ Test operation of group element g on unit vector v

    Here v is a vector in a space, such that v is an 
    instance of class SparseMmVector. g is an element of the 
    group v.space.group operating on the space v.space. 

    The group operation v -> v * g is tested against the operation
    sparse_space(v) * g, where sparse_space(v) is an instance of
    class SparseMmRepVector, which represents the same space as
    v.space, but with vectors given in sparse form.

    If f_mul is set to a function such that f_mul(v, g) = v * g
    holds, an alternative implemetation of the group multiplication
    may be tested. This is useful for some low-level module tests. 
    """
    #print("test_op_case")
    space = MMTestSpace(v.p)
    ref_space = space.ref_space
    v_sp =  ref_space(v)
    if verbose:
        print("p= %d, v=%s, g=%s" % (v.p, v, g))
        print("v_sparse=", v_sp)
    vg = f_mul(v, g) if f_mul else v * g 
    vg.check()
    vg_sp = ref_space(vg) 
    v_sp_g = v_sp * g
    ok = v_sp_g == vg_sp
    if verbose or not ok:
        if not ok:   
            print("p= %d, v=%s, g=%s" % (v.p, v, g))
        if 1:
            print("v * g  expected in sparse format:")
            print(v_sp_g)
            print("v * g  obtained in sparse format:")
            print(vg_sp)            
        print("v*g=", vg) 
        if not ok:   
            raise ValueError("Group operation on vector failed")
    vg_converted = space(vg_sp)
    if vg != vg_converted:
        raise ValueError("Conversion of result from sparse format failed")
 
 

def op_testcases(p):
    space = MMTestSpace(p)
    group = space.group
    for v_tag in "XYZ":  
        for g_tag in "d": 
            for i0 in (1 <<  i for i in range(11)):  
                for i in range(10):
                    g = group(g_tag, 'r')
                    v = space(v_tag, i0, randint(0,23))
                    yield v, g
    for v_tag in "ABCXYZT":  
        for g_tag in "dpxy":   
            for i in range(10):
                v = space(v_tag, 'r')
                g = group(g_tag, 'r')
                yield v, g
        for g_tag in "tl":   
            for e in [1,2,1,2]:
                v = space(v_tag, 'r')
                g = group(g_tag, e)
                yield v, g
 
   
@pytest.mark.mm_op
def test_op(f_mul = None, verbose = 0):
    print("Testing group operation on sparse vectors")
    i = 0
    for p in PRIMES:
        for args in op_testcases(p):
            if verbose: print("Test", i+1)
            one_test_op(*args, f_mul = f_mul, verbose = verbose)
            i += 1
            #print("<%d>" % i, end = "", flush=True)
    print("Test passed")
   
   

########################################################################
# tests with a random vector v
########################################################################


def one_test_rand_op(v, g, basis_vectors, f_mul = None, verbose = 0):
    """Test operation of group element g on random vector v

    We compute v * g with the fast C and with the slow python 
    implementation and check that the results are equal.

    v must be a vector in a space which is an instance of class
    MMTestSpace, and g must be an element of the monster group
    such that v * g can be computed. So g may e.g. be an element
    of v.space.group. If v is an instance of class MMTestSpace
    then a random vector in the space given by v is created.

    Since the python implementation is rather slow, we compute v*g
    in python only for the projection of v to the space spanned 
    by the given 'basis_vectors'. Here a set of basis vectors
    must be encoded in the same way as in function group_blocks().

    If g operates non-monomially on v, it may be necessary to
    consider more basis vectors than given by the argument
    'basis_vectors'. Then function group_blocks() is used to
    extend the set of basis vectors in suitable way if necesary.

    The extension of basis vectors may fail if the structure of
    the matrix of g is too complicated. But that extension is
    always successful if g contains at most one atom or if all
    atoms in g operate monomially. So we can test the operation
    of all atoms in g.

    If f_mul is set to a function such that f_mul(v, g) = v * g
    holds, an alternative implemetation of the group multiplication
    may be tested in te same way as in function test_op_case(). 
    """
    #if isinstance(v, AbstractMmRepSpace):
    #    v = v.rand_vector()
    space = MMTestSpace(v.p)
    ref_space = space.ref_space
    v_g = f_mul(v, g) if f_mul else v * g
    if verbose:
        print("p = %d, g = %s" % (v.p, g))
    basis_block = group_blocks(basis_vectors, g)
    v_sparse = ref_space('S', v.get_sparse(basis_block))
    if verbose:
        print("\nTest operation v * g")
        print("len(v_sparse) =", len(v_sparse))
        print("sparse(v) =", v_sparse)
    v_sparse_g = v_sparse * g
    v_g_sparse = ref_space(v_g.projection(*v_sparse_g.as_tuples()))
    #v_g_sparse = ref_space(v_g.projection(v_sparse_g))
    ok  = v_sparse_g == v_g_sparse
    if verbose or not ok:
        if not ok:
            print("p = %d, g = %s" % (v.p, g))
            print("basis block = ", v.get_sparse(basis_block))
            print("len(sparse v) =", len(v_sparse),
                            ", type =", type(v_sparse))
            print("sparse(v) =", v_sparse)
            print("len(sparse(v*g)) =", len(v_g_sparse))
            print("sparse(v * g) =", v_g_sparse)
        print("len(sparse(v) * g) =", len(v_sparse_g))
        print("sparse(v) * g =", v_sparse_g)
        print("diff =", v_g_sparse - v_sparse_g)
        print("")
    if not ok:
        raise ValueError("Group operation in C differs from python op.")
    
    

@pytest.mark.mm_op
def test_rand_op(n_tests = 3, f_mul = None, verbose = 0):
    print("Testing group operation on random vectors")
    for i in range(n_tests):
        for p in PRIMES:
            space = MMTestSpace(p)
            group = space.group
            for atom in "dpxytl":
                g = group(atom, "n")
                if atom == "l":
                    basis = "D" + "BCTX" * 20 + "A" * 4 + "YZ" * 2
                elif atom == "t":
                    basis = "D" * 3 + "ABC" * 7 +  "T" * 2 + "XYZ" * 20
                else:
                    basis = "D" * 3 + "ABC" * 10 + "TXYZ" * 20 
                for i in range(n_tests): 
                    v = space("R")               
                    one_test_rand_op(v, g, basis, f_mul = f_mul,
                                  verbose = verbose)
    print("Test passed")


########################################################################
# tests operation on tag 'A' of the vector
########################################################################


def tag_A_testwords():
    for atom in "dpxyl":
        yield [(atom, "n")]
    yield [(c, "n") for c in "dpxyl" * 3]


@pytest.mark.mm_op
def test_rand_op_tag_A(n_tests = 4, f_mul = None, verbose = 1):
    print("Testing group operation on random vectors")
    for i in range(n_tests):
        for p in PRIMES:
            space = MMTestSpace(p)
            group = space.group
            for w in tag_A_testwords():
                g = group(w)
                #print(g)
                v = space('R')
                a_g = v.copy()
                v_g = v * g
                len_g = len(g.mmdata)
                res = mm_op_word_tag_A(v.p, a_g.data, g.mmdata, len_g, 1)
                assert res == 0
                #print(v_g['A'] - a_g['A'] + 3)
                assert (v_g['A'] ==  a_g['A']).all()
                res = mm_op_word_tag_A(v.p, v.data, g.mmdata, len_g, -2)
                assert res == 0
                res = mm_op_word_tag_A(v.p, v.data, g.mmdata, len_g, 3)
                assert res == 0
                assert (v['A'] ==  a_g['A']).all()
    print("Test passed")




