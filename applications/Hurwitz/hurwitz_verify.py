import sys
from random import sample, randint, shuffle


sys.path.append(r".")
from find_generators import  MM
from check import mm_product

from hurwitz_monster_samples import hurwitz_monster_data


z = MM('x', 0x1000)

def hurwitz_get_orders(z, g, verifiers):
    z, g = MM(z), MM(g)
    d1 = [z * g, z * g**(-1)]
    d2 = {}
    for i in range(4):
        d2[i] = d1[i & 1] * d1[(i >> 1) & 1]
    orders = []
    for x in verifiers:
        ver = mm_product([d2[(x >> i) & 3] for i in range(0, 64, 2)])
        orders.append(ver.order())
    return orders

def hurwitz_check_sample(z, g, verifiers):
    if z.order() != 2: 
        return False
    if g.order() != 3: 
        return False
    if (z*g).order() != 7: 
        return False
    orders = hurwitz_get_orders(z, g, verifiers)
    if not 94 in orders:
        return False
    for o in [41, 59, 71]:
        if o in orders:
            return True
    return False


def commutator(a, b):
    return (b**-1)**a * b

def compute_commutators():
    s = """Presentation of the Monster as a Hurwitz group with 
generators a, b and relations  a^2 = b^3 = (a*b)^7 = 1\n"""
    print(s)
    s_ok = "Presentation found, commutator [a,b] has order %s"
    i = 0
    for g3, ge, verifiers in hurwitz_monster_data:
        g = MM(g3) ** MM(ge)
        if  hurwitz_check_sample(z, g, verifiers):
            o = commutator(z, g).order()
            print(s_ok %  o)
        i += 1
    print("%d presentations found" % i)


if __name__ == "__main__":
    compute_commutators()
        
 
