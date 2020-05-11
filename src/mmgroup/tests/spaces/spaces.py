from mmgroup.mm_space import MMSpace
from mmgroup.mm_space import characteristics
from mmgroup.mm_group import MMGroup

from mmgroup.tests.spaces.sparse_mm_space import SparseMmSpace
from mmgroup.tests.groups.mgroup_n import MGroupN

spaces_dict = {}

g = MMGroup()
ref_g = MGroupN(kernel = None)
g.set_preimage(ref_g, tuple)
ref_g.set_preimage(g, tuple)

def make_test_space(p):
    sp = MMSpace(p, g)
    ref_sp = SparseMmSpace(p, ref_g)
    sp.ref_group = ref_g
    sp.ref_space = ref_sp
    sp.set_preimage(ref_sp, sp.via_sparse) 
    ref_sp.set_preimage(sp, ref_sp.via_sparse) 
    return sp
       
def MMTestSpace(p):
    global spaces_dict
    try:
        return spaces_dict[p]
    except KeyError:
        spaces_dict[p] = make_test_space(p)
        return spaces_dict[p]


