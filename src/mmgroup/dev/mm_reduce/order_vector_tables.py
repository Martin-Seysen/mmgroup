import os
import sys
import time
import numpy as np
import re


import_pending = True


def import_all():
    global MMV15, import_pending
    global MM0,  MMV, MMVector, MMSpace
    global get_order_vector, order_vector_from_dict
    from mmgroup.structures.mm0_group import MM0
    from mmgroup.mm_space import MMV, MMVector, MMSpace
    from mmgroup.dev.mm_reduce.order_vector import get_order_vector
    from mmgroup.dev.mm_reduce.order_vector import order_vector_from_dict
    MMV15 = MMV(15)
    import_pending = False
  

NAMES = [
   "TAG_DATA", 
   "S_G71", "S_V71", "S_GA", "DIAG_VA", "S_G94", "S_V94",
   "TAGS_Y", "TAGS_X", "TAG_SIGN"
]


def order_vector_from_py_data():
    if import_pending:
        import_all()
    ov, o_tag, data = get_order_vector()
    data[ "TAG_DATA"] = o_tag
    return ov, o_tag, data


def order_vector_from_C_data():
    if import_pending:
        import_all()
    data = {}
    LEN = 128
    from mmgroup.mm_reduce import  mm_order_load_tag_data
    from mmgroup.mm_reduce import  mm_order_load_vector
    for i, name in enumerate(NAMES):
        a = np.zeros(LEN, dtype = np.uint32)
        length = mm_order_load_tag_data(i, a, LEN)
        data[name] = a[:length]
    o_tag = data["TAG_DATA"]
    ov = MMV15(0)
    mm_order_load_vector(ov.data)
    return ov, o_tag, data



def compress_order_vector_data(ov):
    t = []
    v = ov.data
    def compress24(start, length):
        for i in range(start, start + 2*length, 2):
             t.append(int(v[i]) & 0xffffffff)
             t.append(int(v[i]) >> 32)
             t.append(int(v[i+1]) & 0xffffffff)
    def compress64(start, length):
        for i in range(start, start + 4*length):
             t.append(int(v[i]) & 0xffffffff)
             t.append(int(v[i]) >> 32)
    OFS_ABC = 0
    OFS_T = 72*2
    OFS_XYZ = OFS_T + 759*4
    compress24(OFS_ABC, 72)
    compress64(OFS_T, 759)
    compress24(OFS_XYZ, 3*2048)
    return np.array(t, dtype = np.uint32)


class OrderVectorTable:
    directives = {}
    tables = None

    def compute_data(self):
        """Compute data for table classes

        We cannot to this ab initio; because the required inputs 
        are not available if sys.arg contains tan argument ``mockup``.
        """
        if not self.__class__.tables is None:
            return
        if import_pending:
            import_all()
        ov, tag, data =  order_vector_from_py_data()  
        a_ov =  compress_order_vector_data(ov)
        tag_data = []
        tag_indices = [0]
        for i, name in enumerate(NAMES):
            a0 = list(data[name])
            tag_data += a0
            tag_indices.append(tag_indices[-1] + len(a0))
        tag_data = np.array(tag_data, dtype = np.uint32)
        tag_indices = np.array(tag_indices, dtype = np.uint16)
        self.__class__.tables = {
            "ORDER_VECTOR_DATA": a_ov,
            "ORDER_VECTOR_TAG_DATA": tag_data,
            "ORDER_VECTOR_TAG_INDICES": tag_indices,
            "ORDER_VECTOR_NUM_TAG_DATA": len(NAMES),
            }
    
    def __init__(self, *args):
        self.compute_data()


class Mockup_OrderVectorTable:
    a_ov =  np.array([0], dtype = np.uint32)
    tag_data = np.array([0], dtype = np.uint32)
    tag_indices = np.array([0,0], dtype = np.uint16)
    tables = {
        "ORDER_VECTOR_DATA": a_ov,
        "ORDER_VECTOR_TAG_DATA": tag_data,
        "ORDER_VECTOR_TAG_INDICES": tag_indices,
        "ORDER_VECTOR_NUM_TAG_DATA": 1
    }
    directives = {}

    def __init__(self, *args):
        pass



def check_order_vector_data():
    if import_pending:
        import_all()
    from mmgroup.mm_reduce import  mm_order_compare_vector
    from random import randint

    ov, o_tag, data = order_vector_from_C_data()
    ov_ref, o_tag_ref, _ = order_vector_from_dict(data)
    ov, ov_ref = ov.data, ov_ref.data
    if not np.array_equal(ov, ov_ref):
        for i, w in enumerate(ov):
            w_ref = ov_ref[i]
            if w != w_ref:
                print("Order vector error at index %d" % i)
                print("Found: %s, expected: %s" % (hex(w), hex(w_ref)))
                raise ValueError("Error in order vector")
    assert np.array_equal(o_tag, o_tag_ref)
    if (mm_order_compare_vector(ov_ref)):
       err = "Comparing with order vector failed" 
       raise ValueError(err) 
    positions = [(0,0), (72*2, 63), (72*2 + 759*4 - 1, 35)]
    for i in range(30, len(ov) - 1, 537):
        bmax = 63 if 72*2 < i < 72*2 + 759*4 else 31
        positions.append( (i, randint(0,bmax)) ) 
    for i, sh in positions:
         a = np.copy(ov_ref)
         a[i] = int(a[i]) ^ (1 << sh)
         result = mm_order_compare_vector(a)  
         if (result != 1):
             err = "Comparing of bad vector with order vector failed" 
             raise ValueError(err) 
              


