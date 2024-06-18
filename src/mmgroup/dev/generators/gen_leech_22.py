"""Compute tables für functions in file gen_reduce_22.c"""


import sys
import random
if __name__ == "__main__":
    sys.path.append("../../..")

from mmgroup.dev.mat24.mat24_ref import Mat24 as m



def iter_octads_22():
    """Iterate over the 77 octads containing entries 2 and 3

    The function yields the 77 octads o containing the entries
    2 and 3 as bit vectors. It also yields a list bl of five
    entries of octad o, excluding the entries 2 and 3.
    So the function yields 77 pairs (o, bl).

    The computed pairs satisfy the following additional conditions.
    Precisely the first 21 octads contain an entry 1. None of the
    computed bit lists bl contains an entry 1.
    """
    n = k = 0
    tail = []
    for i in range(759):
        o = int(m.octad_to_vect(i))
        if o & 0xc == 0xc:
            n += 1
            length, bl = m.vect_to_bit_list(o & ~0xc)
            assert length == 6
            bl = bl[:length]
            if 1 in bl:
                bl.remove(1)
                yield o, bl
                k += 1
            else:
                tail.append((o, bl[:5]))
    assert n == 77, n
    assert k == 21, k
    yield from tail



def iter_vectors_22():
    """Yield Leech lattice mod 2 vectors from function iter_octads_22

    For each of the 77 octads obtained from function iter_octads_22()
    there are 16 type-2 vectors v2 in the Leech lattice mod 2 such
    v2 + w2  is also of type 2. Here w2 is the vector
    (0, 0, 4, -4, 0,...,0) in the Leech lattice mod 2.
    The function returns quintuples (x1, x2, x3, x4, ov)
    such that  v = ov + u1*x1 + u2*x2 + u3*x3 + u4*x4 + u5*Omega
    is a type_2 vector as described above for any u1, u2, u3, u4
    in [0, 1]. Here u5 depends on  u1, u2, u3, u4 as follows:
    u5 is 1 if u1 + u2 + u3 + u4 in [0, 1, 4]; and u5 = 0
    otherwise. Here x1,...,x4 are even cocode vectors.

    Precisely the vectors v computed from the first 21 octads satisfy
    the following additional conditions. v2 + w1  is also of type 2;
    and  v2 + w1 + w2 is of type 4, for w1  = (0, 4, -4, 0,...,0);
    """
    for o, bl in iter_octads_22():
        og = m.vect_to_gcode(o)
        ov = (og << 12) + (int(m.theta_table[og & 0x7ff]) & 0xfff)
        ov ^= m.vect_to_cocode((1 << 2) | (1 << bl[0]))
        lst = []
        for i in range(1, 5):
            lst.append(m.vect_to_cocode((1 << bl[0]) | (1 << bl[i])))
        yield lst + [int(ov)]

def make_table():
    tbl_gcode, tbl_cocode = [], []
    for lst in iter_vectors_22():
       tbl_gcode += lst[4:]
       tbl_cocode += lst[:4]
    return tbl_gcode, tbl_cocode





def w_mod4():
    """Compute a certain parity function

    The function returns a value a such that bit i of a is
    set if and only if the integer i has bit parity 0, 1, or 4;
    for 0 <= i < 16.

    So this function may be used to compute the bit u5 defined
    in function iter_vectors_22().
    """
    a = 0
    for i in range(16):
        if m.bw24(i) in [0,1,4]:
            a |= 1 << i
    return a



def make_odd_coc_table():
    r"""Compute a certain table of odd  covectors

    the function returns a table of 22 covectors x_delta with
    delta = {i, 2, 3}; 0 <= i < 24; i != 2, 3.
    """
    coc_table = []
    for i in list(range(4, 24)) + [0, 1]:
        vector = (1 << i) ^ 0xc
        coc_table.append(m.vect_to_cocode(vector))
    return coc_table



class Prime4600:
    primes = [p for p in range(5, 75, 2) if p % 5 and p % 7]
    p = m = 0

    @classmethod
    def compute_p(cls):
        start = (4600 // 12) * 12 + 11
        for p in range(start, start+800, 12):
            ok = 1
            for pi in cls.primes:
                ok = p % pi
                if not ok:
                    break
            if ok:
                cls.p = p
                return p
        raise ValueError("No suitable prime found")

    @classmethod
    def factors(cls, n):
        f = []
        for p in [2,3] + cls.primes:
            while n % p == 0:
                f.append(p)
                n = n // p
        if n > 1:
            f.append(n)
        return f

    @classmethod
    def cofactors(cls, n):
        return [n // x for x in set(cls.factors(n))]

    @classmethod
    def multiplier(cls):
        if not cls.p:
            cls.compute_p()
        p = cls.p
        cofactors = cls.cofactors(p - 1)
        phi = (5 ** 0.5 - 1) / 2
        start = int(phi*cls.p)
        def iter_dist():
            yield start
            for i in range(1, 100):
                yield start + i
                yield start - i
        for m in iter_dist():
            assert pow(m, p-1, p) == 1
            ok = True
            for e in cofactors:
                ok = pow(m, e, p) != 1
                if not ok:
                    break;
            cls.m = m
            return m
        raise ValueError("No suitable multiplier found")

    @classmethod
    def check(cls, verbose = 1):
        if not cls.m:
            cls.multiplier()
        p, m = cls.p, cls.m
        s = set()
        x = 1
        for i in range(p-1):
            x = x * m % p
            s.add(x)
        assert s == set(range(1,p))
        if verbose:
            print("Recommended prime: %d, multplier: %d" % (p, m))


class Tables:
    w = w_mod4()
    odd_coc_table = make_odd_coc_table()
    tbl_gcode, tbl_cocode = make_table()


    tables = {
       "GenLeech_v22_table_gcode" : tbl_gcode,
       "GenLeech_v22_table_cocode" : tbl_cocode,
       "GenLeech_v22_weights" : w,
       "GenLeech_v22_odd_coc_table" : odd_coc_table,
    }

    directives = {}


if __name__ == "__main__":
    Prime4600.check()

