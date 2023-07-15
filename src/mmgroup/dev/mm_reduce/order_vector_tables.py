import os
import sys
import time
import numpy as np
import re

if __name__ == "__main__":
   sys.path.append(os.path.join('..', '..', '..'))


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


#####################################################################
# Relations in Frobenius group of order 21 generated by tags 't', 'l'
#####################################################################

def tlt_element(n, tags):
    a = [(t, 1 + ((n >> j) & 1)) for j, t in enumerate(tags)]
    return MM0(a)


def relations_tlt(tags, v):
    d = {}
    for n in range(1 << len(tags)):
        d[n] = (v * tlt_element(n, tags)).hash()
    return d

def tlt_to_ltl(verbose =  0):
    if verbose:
        s = "Make table with map  t^i l^j t^k -> l^i' t^j' l^k'" 
        print(s)
    if import_pending:
        import_all()
    for _ in range(100):
        v = MMV15('R')
        src = relations_tlt("tlt", v)
        dest = relations_tlt("ltl", v)
        h_src, h_dest = src.values(), dest.values()
        if len(h_src) == 8 and set(h_src) == set(h_dest):
            rev_dest = dict([(h, k) for k, h in dest.items()])
            d_out = dict([(k, rev_dest[src[k]]) for k in src]) 
            v1 = MMV15('R')
            for j in range(8):
                tlt = tlt_element(j, "tlt")
                ltl = tlt_element(d_out[j], "ltl")
                assert v1* tlt == v1 *ltl
            if verbose:
                print(d_out)
            return d_out
        if verbose:
            print("Attempt failed")
    raise ValueError("tlt->ltl tag conversion has failed")

def tlt_conversion_table():
    d = tlt_to_ltl()
    return sum(d[i] << (4*i) for i in range(8))

#####################################################################
# End of dealing with in Frobenius group of order 21
#####################################################################



class OrderVectorTable:
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
        ov, tag, data =  order_vector_from_py_data()
        ov_hash = ov.hash()  
        a_ov =  compress_order_vector_data(ov)
        tag_data = []
        tag_indices = [0]
        for i, name in enumerate(NAMES):
            a0 = list(data[name])
            tag_data += a0
            tag_indices.append(tag_indices[-1] + len(a0))
        tag_data = np.array(tag_data, dtype = np.uint32)
        tag_indices = np.array(tag_indices, dtype = np.uint16)
        tlt_conversion = tlt_conversion_table()
        self.__class__.tables = {
            "ORDER_VECTOR_DATA": a_ov,
            "ORDER_VECTOR_TAG_DATA": tag_data,
            "ORDER_VECTOR_TAG_INDICES": tag_indices,
            "ORDER_VECTOR_NUM_TAG_DATA": len(NAMES),
            "ORDER_VECTOR_HASH": ov_hash,
            "ORDER_VECTOR_TLT_CONVERSION": tlt_conversion,
            }
    
    def __init__(self, *args, **kwds):
        self.compute_data()


class Mockup_OrderVectorTable:
    def __init__(self, *args, **kwds):
        pass
    a_ov =  np.array([0], dtype = np.uint32)
    tag_data = np.array([0], dtype = np.uint32)
    tag_indices = np.array([0,0], dtype = np.uint16)
    tables = {
        "ORDER_VECTOR_DATA": a_ov,
        "ORDER_VECTOR_TAG_DATA": tag_data,
        "ORDER_VECTOR_TAG_INDICES": tag_indices,
        "ORDER_VECTOR_NUM_TAG_DATA": 1,
        "ORDER_VECTOR_HASH": 0,
        "ORDER_VECTOR_TLT_CONVERSION": 0,
    }
    directives = {}


    


Tables = OrderVectorTable
MockupTables = Mockup_OrderVectorTable



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
              

if __name__ == "__main__":
   tlt_to_ltl(verbose = 1)
   print ("Hex table:", hex(tlt_conversion_table()))

