"""Compute tables für functions in file gen_reduce_22.c"""

from mmgroup.dev.mat24.mat24_ref import Mat24 as m



def iter_octads_22():
    """Iterate over the 77 octads containing entries 2 and 3

    The function yields the 77 octads o containing the entries
    2 and 3 as bit vectors. It also yields a list bl of five
    entries of octad o, excludin the entries 2 and 3.
    So the function yields 77 peits (o, bl)
    """
    n = 0
    for i in range(759):
        o = m.octad_to_vect(i) 
        if o & 0xc == 0xc:
            n += 1
            bl = m.vect_to_bit_list(o & ~0xc)
            #print(hex(o), bl[1][:5])
            yield int(o), bl[1][:5]
    assert n == 77, n



def iter_vectors_22():
    """Yield Leech lattice mod 2 vectors from function iter_octads_22

    For each of the 77 octads obtained from function iter_octads_22()
    there are 16 type-2 vectors v2 in the Leech lattice mod 2 such
    v2 + w2  is also of type 2. Here w2 is the vector
    (0, 0, 4, -4, 0,...,0) in the Leech lattice mod 2.
    The function returns quintuples (x1, x2, x3, x4, ov)
    such that ov + u1*x1 + u2*x2 + u3*x3 + u4*x4 + u5*Omega is
    a type_2 vector as described above for any u1, u2, u3, u4
    in [0, 1]. Here u5 depends on  u1, u2, u3, u4 as follows:
    u5 is 1 if u1 + u2 + u3 + u4 in [0, 1, 4]; and u5 = 0
    otherwise. Here x1,...,x4 are even cocode vectors.
    """
    for o, bl in iter_octads_22():
        og = m.vect_to_gcode(o)
        ov = (og << 12) + (m.theta_table[og & 0x7ff] & 0xfff)
        ov ^= m.vect_to_cocode((1 << 2) | (1 << bl[0]))
        lst = []
        for i in range(1, 5):
            lst.append(m.vect_to_cocode((1 << bl[0]) | (1 << bl[i])))
        #print([hex(x) for x in lst + [ov]])
        yield lst + [int(ov)]


def compress(x):
    """Compress element obtained from iter_vectors_22()

    A 24-bit integer obtained from function iter_vectors_22()
    decribes an element x of the Leech lattice mod 2. By
    construction, bits 9 and 11 of such a value x are always zero.
    We compress x by mapping bits (23,...,12,10,8,...,0) of x
    to bits (21,...,10,9,8,...,0). We return that compresses value.
    """
    return (x & 0x1ff) + ((x & 0x400) >> 1) + ((x & 0xfff000) >> 2)


def make_table():
    """Make table from function iter_vectors_22()

    The function creates a table of 77 64-bit integers. Any such
    integers correspond to a (compressed) quintuple
    (x1, x2, x3, x4, ov) as obtained by function iter_vectors_22().
    Here function compress() is used to compress x1, x2, x3, x4
    from 12 to 10 bit and to compress ov from 24 to 22 bits.

    The value x_i is stored in bits 10*i+9,...,10*i from an entry;
    and ov is stored in bits 61,...,40 of an entry.
    """
    tbl = []
    for lst in iter_vectors_22():
        a = 0
        for i, v in enumerate(lst):
            #print(i, hex(v), hex(compress(v)))
            a = a + (compress(v) << (10*i))
            #print(".",hex(a))
        tbl.append(a)
        #print(hex(a))
    return tbl


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

class Tables:
    tbl = make_table()
    w = w_mod4()
    #print(tbl)
    #print(hex(w))

    tables = {
       "GenLeech_v22_table" : tbl,
       "GenLeech_v22_weights" : w,
    }

    directives = {}





