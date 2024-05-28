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
    return v.last_timing / iterations


def bench_nobreak(p, operation = [], iterations = ITERATIONS_NOBREAK):
    return bench(p, operation, iterations, break_g = False) 

def bench_weights(frequencies, timings):
    t = 0.0
    for tag, f in frequencies.items():
        t += f * timings[tag] 
    return t   


def quot_ms(f, *args):
    t = f(*args)
    return t, "%9.6f" % (1000 * t)


# frequencies of operators in random mmgroup element
FREQ = {'xy':1, 'p':16, 'l':17, 't':6.25}

MM_REDUCE_OPS = {3:3.0, 15:(3.0 + 2 + 3.0/7.0)}

FREQ_G = {'xy':1, 'p':4, 'l':3}

FREQ_N = {'xy':1, 'p':1}


def benchmark(few = True):
    print("""
Benchmarking monster operations on a %d-bit system.
All times are given in milliseconds.
""" % (INT_BITS)
    )
    p_values = [3, 15] if few else characteristics()
    runtimes = {}
    runtimes_mm = {}
    for p in p_values:
        print("Operation modulo", p)
        op = [('p', 22), ('d', 127)]
        print ("p    ", quot_ms(bench_nobreak, p, op)[1], " local optimization")
        print ("p    ", quot_ms(bench, p, op)[1])

        op = [('p', 23), ('d', 12745645)]        
        runtimes['p'], msg = quot_ms(bench, p, op)
        print ("p odd", msg)

        op = [('x', 1237), ('y', 567),]
        runtimes['xy'], msg = quot_ms(bench, p, op)
        print ("xy   ", msg)

        op = [('l', 2)]
        runtimes['l'], msg = quot_ms(bench, p, op)
        print ("l    ", msg)

        op = [('t', 2)]
        runtimes['t'], msg = quot_ms(bench, p, op)
        print ("t    ", msg)

        _, msg = quot_ms(bench_weights, FREQ_N, runtimes)
        print ("N_x0 ", msg)

        _, msg = quot_ms(bench_weights, FREQ_G, runtimes)
        print ("G_x0 ", msg)

        runtimes_mm[p], msg = quot_ms(bench_weights, FREQ, runtimes)
        print ("MM   ", msg)

        if few:
           continue
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

    _, t_mm = quot_ms(bench_weights, MM_REDUCE_OPS, runtimes_mm)
    print("Estimated time for group operation:", t_mm)

@pytest.mark.mm_op
@pytest.mark.slow
@pytest.mark.bench
@pytest.mark.user
def test_benchmark():
    benchmark()


if __name__ == "__main__":
    benchmark()


