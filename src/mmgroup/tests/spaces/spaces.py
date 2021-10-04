from mmgroup.mm_space import MMSpace, MMV, MMVector
from mmgroup.mm_space import characteristics
from mmgroup.structures.mm0_group import MM0Group

from mmgroup.tests.spaces.sparse_mm_space import SparseMmSpace
from mmgroup.tests.spaces.sparse_mm_space import SparseMmV
from mmgroup.tests.groups.mgroup_n import MGroupN


#print("module mmgroup.tests.spaces.spaces is deprecated!!")


spaces_dict = {}

g = MM0Group()
ref_g = MGroupN()


class TestSpace:
    def __init__(self, p):
        self.p = p
        self.space = MMV(p) 
        self.ref_space = SparseMmV(p)
        self.group = g     
        self.ref_group = ref_g 
    def __call__(self, *args):
        return self.space(*args)



def MMTestSpace(p):
    global spaces_dict
    try:
        return spaces_dict[p]
    except KeyError:
        spaces_dict[p] = TestSpace(p)
        return spaces_dict[p]


