r"""This little program is yet to be documented.

We consider the automorphism ``s`` of the Y_555 diagram that exchanges
b_2 with b_3, c_2 with c_3, d_2 with d_3, etc. We will show that the
image ``img_s`` of ``s`` in the Bimonster is an element of the diagonal 
of the the direct product of the Monster with itself. Here the two
factors are 2A involutions in the Monster.

"""


from collections import defaultdict
from random import sample
from multiprocessing import Pool

try:
    import mmgroup
except (ImportError, ModuleNotFoundError):
    # get the mmgroup package from its inplace location it not installed
    import os
    import sys
    sys.path.append(os.path.join('..', '..', 'src'))
    import mmgroup


from mmgroup.bimm import BiMM, AutP3, P3_BiMM, AutP3_BiMM, P3_node
from mmgroup import MM

def legswap_involution_and_generators(verbose =  0):
    # Automorphism ``s`` is uniquely defined by imagess of a, c1, c2, and c3
    s = AutP3('a:a, c1:c1, c2:c3, c3:c2')
    text = """Let 's' be the automorphism of the Y_555 graph that exchanges
b_2 with b_3, c_2 with c_3, d_2 with d_3, etc., and fixes a, b1, c1, etc.
Let 'img_s' be the image of 's' in the Bimonster.
"""
    if verbose: print(text)
    # Let ``img_s`` be the image of ``s`` in he Bimonster
    img_s = AutP3_BiMM(s)
    if verbose: print("'img_s' is ", img_s)
    # Decompose the image ``img_s``
    m1, m2, e = img_s.decompose()
    # Check that 'img_s' is in the direct square of the Monster 
    assert e == 0
    # Check that ``img_s`` is in the diagonal 
    assert m1 == m2
    if verbose: print("'img_s' is in the diagonal of the direct square of the Monster.")
    I1, _ =  m1.conjugate_involution()
    assert I1 in [1,2]
    class_ = '2A' if I1 == 1 else '2B'
    msg = "The factors of 'img_s' are %s involutions in the Monster."
    if verbose: print(msg % class_) 
    gen  = [['a']]
    for x in "bcde":
        gen.append(['%s1' % x])
        gen.append(['%s2' % x, '%s3' % x])
    for y in gen + [['f1'],['f2', 'f3']]:
        #print(y)
        assert P3_BiMM(y) == P3_BiMM(y)**img_s
    return s, gen
    

def sample_BiMM(generators, length = 400):
    """Sample a word of (the even part of) a subgroup of the Bimonster

    The function returns a random word of the even part G of the 
    subgroup of the Bimonster generated by the Coxeter relections
    
    a, b1, c1, d1, e1, b2*b3, c2*c3, d2*d3, e2*e2.
    """ 
    while True:
        g = sample(generators*length, length)
        g_flat = [item for sublist in g for item in sublist]
        if len(g_flat) & 1 == 0:
            return P3_BiMM(g_flat)


    
global s, gen
s = None
gen = None
s, gen = legswap_involution_and_generators()

def find_orders(nsamples):
    """Compute orders of element of a subgroup of the Bimonster

    The function samples ``nsamples`` elements of the subgroup G of 
    the Bimonster described in function ``sample_BiMM``. All these 
    elements are in the direct product of the Monster with itself.
    
    From each sample ``m`` we take its two components ``m1, m2`` 
    in the Monster. Let ``s`` be the leg-swapping automorphism 
    as in function ``legswap_involution_and_generators``, and 
    let ``s1, s2`` the decomposition of ``s`` in the direct product
    of the Monster with itself. Since each sample ``m`` in G  
    commutes with ``s``, its components ``m1`` and ``m2`` commute 
    with ``s1`` and ``s2``, respectively. Note that ``s1`` and 
    ``s2`` are 2A involutions in the Monster. Thus ``m1`` and ``m2``
    are in the centralizers ``C1`` and ``C2`` of ``s1`` and 
    ``s2``, resspectively. Both, ``C1`` and ``C2`, have structure 
    2.B, where B is the Baby Monster.

    The function computes the orders ``o1`` and ``o2`` of 
    ``m1`` and ``m2`` in the factor groups ``C1/2`` and ``C2/2`` 
    obtained from ``C1`` and ``C2`` by factoring out their 
    centers of order 2. Clearly, ``C1/2`` and ``C2/2`` are
    isomorphic to the Baby Monster. 

    We conjecture that G is the direct product of 2.B with itself.
    So by observing the distributions of the orders ``o1`` and 
    ``o2`` we may see that the projection of G to any of the 
    factors of the direct product of the Monster with itself
    is 2.B.
         
    The function returns a dictionary. The keys are the pairs 
    ``(o1, o2)`` of orders found. The value of a key is the 
    number of elements of G found for that pair of orders. 
    """
    s1, s2, _ = AutP3_BiMM(s).decompose()
    orders = defaultdict(int)
    for i in range(nsamples):
        g = sample_BiMM(gen)
        m1, m2, e = g.decompose()
        assert e == 0
        order1, h1 = m1.half_order() 
        if h1 == s1:
            order1 = order1 // 2
        order2, h2 = m2.half_order() 
        if h2 == s2:
            order2 = order2 // 2
        #print(i, nsamples, order1, order2)
        orders[(order1, order2)] += 1
    return orders



def display_orders(orders, index):
    txt = "Distribution of orders of component %d of the group"
    print(txt % index)
    d = defaultdict(int)
    for o, n in orders.items():
        d[o[index]] += n
    for o in sorted(d):
        print("%s: %3d" %  (o, d[o]))
    
    
def display_order_pairs(orders):
    t = """
Orders of samples of elements of component 1 of (the even part of)
the group generated by a, b1, c1, d1, e1, b2*b3, c2*c3, d2*d3, e2*e2.
We display a few selected pairs or orders only."""
    print(t)
    for o in sorted(orders):
        disp = not 0.5 <= o[1]/o[0] <= 2
        disp |= o[0] == 47 and o[1] != o[0]
        disp |= o[1] == 47 and o[1] != o[0]
        if disp:
            print("%s: %3d" %  (o, orders[o]))




if __name__ == "__main__":
    legswap_involution_and_generators(verbose =  1)
    with Pool() as pool:
        results = pool.map(find_orders, [100] * 10)
        #results = pool.map(find_orders, [10] * 4)
        pass
    pool.join()
    orders = defaultdict(int)
    for d in results:
        for o, n in d.items():
            orders[o] += n
    orders = dict(orders)
    display_order_pairs(orders)
    for index in range(2):
        display_orders(orders, index)
    


