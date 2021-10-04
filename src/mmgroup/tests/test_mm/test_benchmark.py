import sys
import os

sys.path.insert(0, os.path.abspath("src"))
import mmgroup

import pytest


from mmgroup import MMSpace, characteristics, INT_BITS, MMV, MM0

#space = dict()
#for p in characteristics():
#    space[p] = MMSpace(p)
group = MM0

ITERATIONS = 10000
ITERATIONS_NOBREAK = 1000000



def bench(p, operation = [], iterations = ITERATIONS, break_g = True):
    v = MMV(p)('R')
    #print(operation)
    g = MM0(operation)
    v.mul_exp(g, iterations, break_g = break_g)
    return v.last_timing, iterations


def bench_nobreak(p, operation = [], iterations = ITERATIONS_NOBREAK):
    return bench(p, operation, iterations, break_g = False) 

def quot(f, *args):
    q, r = f(*args)
    return q/r

def quot_ms(f, *args):
    return "%9.6f" % (1000 * quot(f, *args))





def benchmark():
    print("""
Benchmarking monster operations a %d-bit system.
All times are given in milliseconds.
""" % (INT_BITS)
    )    
    for p in characteristics():    
        print("Operation modulo", p)
        op = [('p', 22), ('d', 127)]
        print ("p    ", quot_ms(bench_nobreak, p, op), " local optimization")

        print ("p    ", quot_ms(bench, p, op))
        op = [('p', 23), ('d', 12745645)]
        print ("p odd", quot_ms(bench, p, op))

        op = [('x', 1237), ('y', 567),]
        print ("xy   ", quot_ms(bench, p, op))

        op = [('l', 2)]
        print ("l    ", quot_ms(bench, p, op))

        op = [('t', 2)]
        print ("t    ", quot_ms(bench, p, op))

        op = [('t', 2), ('p', 23), ('d', 12745645),
            ('x', 1237), ('y', 567),
            ('l', 1),  ('p', 442),
            ('l', 1),  ('p', 142345345),
            ('l', 1),  ('p', 44865756), ('y', 567),
        ]

        print ("elem ", quot_ms(bench, p, op))


        op = [('t', 2), ('p', 23), ('d',12745645),
            ('x', 1237), ('y', 567),
            ('l', 1)
        ]
        print ("rand ", quot_ms(bench, p, op))
        print("")


@pytest.mark.very_slow
@pytest.mark.slow
@pytest.mark.bench
@pytest.mark.user
def test_benchmark():
    benchmark()


if __name__ == "__main__":
    benchmark()


