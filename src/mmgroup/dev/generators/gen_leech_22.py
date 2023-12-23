from mmgroup.dev.mat24.mat24_ref import Mat24 as m



def iter_octads_22():
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
    return (x & 0x1ff) + ((x & 0x400) >> 1) + ((x & 0xfff000) >> 2)


def make_table():
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





