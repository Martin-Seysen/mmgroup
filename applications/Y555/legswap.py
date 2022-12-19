r"""This little program is yet to be documented.

We consider the automorphism ``s`` of the Y_555 diagram that exchanges
b_2 with b_3, c_2 with c_3, d_2 with d_3, etc. We will show that the
image ``img_s`` of ``s`` in the BiMonster is an element of the diagonal 
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
    # Automorphism ``s`` is uniquely defiend by imagess of a, c1, c2, and c3
    s = AutP3('a:a, c1:c1, c2:c3, c3:c2')
    text = """Let 's' be the automorphism of the Y_555 graph that exchanges
b_2 with b_3, c_2 with c_3, d_2 with d_3, etc., and fixes a, b1, c1, etc.
Let 'img_s' be the image of 's' in the BiMonster.
"""
    if verbose: print(text)
    # Let ``img_s`` be the image of ``s`` in he BiMonster
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
    for x in "bcdef":
        gen.append(['%s1' % x])
        gen.append(['%s2' % x, '%s3' % x])
    for y in gen:
        #print(y)
        assert P3_BiMM(y) == P3_BiMM(y)**img_s
    return s, gen
    

def sample_BiMM(generators, length = 400):
    """Sample a word of (the even part of) a subgroup of the BiMonster

    The function returns a random word of (the even part of) the 
    subgroup of the Bimonster generated by the Coxeter relections
    
    a, b1, c1, d1, e1, f1, b2*b3, c2*c3, d2*d3, e2*e2, f2*f3.
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
    """Compute orders of element of a subgroup of the Monster

    The samples ``nsamples`` elements of the subgroup of the
    BiMonster described in function ``sample_BiMM``. All these 
    elements are in the direct product of the Monster with itself.
    The function computes the order of the first component of
    each element. 

    It returns a dictionary. The keys are the computes orders of
    the elements and the values are the nubers of elements of
    a given order that have been found. 
    """
    s1, s2, _ = AutP3_BiMM(s).decompose()
    orders = defaultdict(int)
    for i in range(100):
        g = sample_BiMM(gen)
        m1, m2, e = g.decompose()
        assert e == 0
        order, h = m1.half_order() 
        if h == s1:
            order = order // 2
        orders[order] += 1
    return orders



if __name__ == "__main__":
    legswap_involution_and_generators(verbose =  1)
    with Pool() as pool:
        results = pool.map(find_orders, [100] * 10)
    pool.join()
    orders = defaultdict(int)
    for d in results:
        for o, n in d.items():
            orders[o] += n
    orders = dict(orders)
    text = """Orders of samples of elements of component 1 of (the even part
of) the group generated by a, b1, c1, d1, e1, f1, b2*b3, c2*c3, 
d2*d3, e2*e2, f2*f3."""
    print(text)
    for o in sorted(orders):
        print("%2d: %3d" %  (o, orders[o]))
    


