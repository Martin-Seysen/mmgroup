"""The module creates table for finding the LSN and MSB for an uint64_t

"""

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

from random import randint


if __name__ == "__main__":
    import sys
    import os
    sys.path.append(os.path.join('..','..','..')) 
    from mmgroup.bitfunctions import bitlen, v2
    #from mmgroup.generate_c import UserDirective


#######################################################################
# Find suitable multiplier for obtaining the bit length
#######################################################################


def to_high_bit_64(x):
    x &= 0xffffffffffffffff
    x |= x >> 32
    x |= x >> 16
    x |= x >> 8
    x |= x >> 4
    x |= x >> 2
    x |= x >> 1
    return x


def len_to_high_bit_64(bitlen):
    #print(bitlen)
    return to_high_bit_64(1 << (bitlen-1)) if bitlen else 0


def make_hi_bit64_table(mult):
    hi_table = [0xff] * 128
    for i in range(65):
        x = (1 << i) - 1      
        ind = ((x * mult) >> 57) & 0x7f
        if hi_table[ind] != 0xff:
           #print(hex(ind), i,  table[ind], hex(mult), hex(x))
           return None
        hi_table[ind] = i
    return hi_table


def make_lo_bit64_table(mult):
    lo_table = [0xff] * 128
    for i in range(65):
        x = (1 << i) & 0xffffffffffffffff
        ind = ((x * mult) >> 57) & 0x7f
        if lo_table[ind] != 0xff:
           return None
        lo_table[ind] = i
    return lo_table


def find_hi_bit64_tables():
    for i in range(1, 10000000):
        mult = randint(1, (1 << 64) - 1) | 1
        table = make_hi_bit64_table(mult)
        if table: 
            print("\rHigh multiplier", hex(mult), ",", i, "trials")
            return mult, table 
        else:
            if i & 0xfff == 0: print('.', end = "", flush = True)


def find_lo_bit64_tables():
    for i in range(1, 10000000):
        mult = randint(1, (1 << 64) - 1) | 1
        table = make_lo_bit64_table(mult)
        if table: 
            print("\rLow multiplier", hex(mult), ",", i, "trials", " "*40)
            return mult, table 
        else:
            if i & 0xfff == 0: print('.', end = "", flush = True)




# Results of functions find_hi_bit64_tables and find_hi_bit64_tables:
HIGH_MULTIPLIER = 0xb7c2ad8bd12cd265 #  62262 trials
LOW_MULTIPLIER =  0x12e91e16a99fdf2b #   2962 trials




#######################################################################
# Class Bit64Tables
#######################################################################


class Bit64Tables(object):
    def __init__(self, *args, **kwds):
       pass

    HIGH_MULTIPLIER = HIGH_MULTIPLIER
    LOW_MULTIPLIER = LOW_MULTIPLIER
    HIGH_TABLE = make_hi_bit64_table(HIGH_MULTIPLIER)
    LOW_TABLE = make_lo_bit64_table(LOW_MULTIPLIER)

 
    tables = {
           "BIT64_HIGH_MULTIPLIER" : HIGH_MULTIPLIER,
           "BIT64_LOW_MULTIPLIER" : LOW_MULTIPLIER,
           "BIT64_HIGH_TABLE" : HIGH_TABLE,
           "BIT64_LOW_TABLE" : LOW_TABLE,
    }

    directives = {}


    @classmethod
    def uint64_bitlen(cls, x):
        assert 0 <= x < 0x10000000000000000
        x |= x >> 32;
        x |= x >> 16;
        x |= x >> 8;
        x |= x >> 4;
        x |= x >> 2;
        x |= x >> 1;
        x = x * cls.HIGH_MULTIPLIER;
        return cls.HIGH_TABLE[(x >> 57) & 0x7f]
    
    @classmethod
    def uint64_lowest_bit(cls, x):
        assert 0 <= x < 0x10000000000000000
        x &= -x;
        x = x * cls.LOW_MULTIPLIER;
        return cls.LOW_TABLE[(x >> 57) & 0x7f]


Tables = Bit64Tables


#######################################################################
# Self test
#######################################################################


MASK = 0xffffffffffffffff

def get_testcases():
    for i in range(10):
        yield i
        yield MASK - i
    for l in range(4, 65):
        for d in (0, 1, 2, 3):
            yield ((1 << l) - d) & MASK
        for _ in range(20):
            yield randint(1 << (l-1), (1 << l) - 1)
    for i in range(0, 64):
        for _ in range(20):
            yield ((randint(1, MASK) | 1) << i) & MASK
           

def do_test_high(verbose = 0):
    if verbose: print("\nTesting method uint64_bitlen")
    f = Bit64Tables.uint64_bitlen
    for x in get_testcases():
        bl0 = bitlen(x)
        bl1 = f(x)
        assert bl0 == bl1, (hex(x), bl0, bl1)
        if verbose: print("%20s %2d %2d" % (hex(x), bl0, bl1))


def do_test_low(verbose = 0):
    if verbose: print("\nTesting method uint64_lowest_bit")
    f = Bit64Tables.uint64_lowest_bit
    for x in get_testcases():
        bl0 = v2(x) if x else 64
        bl1 = f(x)
        assert bl0 == bl1, (hex(x), bl0, bl1)
        if verbose: print("%20s %2d %2d" % (hex(x), bl0, bl1))


def test_Bit64Tables(verbose = 0):
    print("Testing class Bit64Tables")
    do_test_high(verbose)
    do_test_low(verbose)
    print("passed")


if __name__ == "__main__":
    test_Bit64Tables(verbose = 1)

