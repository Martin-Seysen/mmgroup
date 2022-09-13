import sys
import os
import time
from datetime import datetime
from collections import OrderedDict, defaultdict
from random import sample
from multiprocessing import Pool

import numpy as np


try:
    import mmgroup
    from mmgroup import Xsp2_Co1, MMV, MM0, MM
except (ImportError, ModuleNotFoundError):
    sys.path.append(os.path.join('..', '..', 'src'))
    import mmgroup
    from mmgroup import Xsp2_Co1, MMV, MM0, MM
  
sys.path.append(r".")


########################################################################
# Number of CPUs used in parallel. NPROCESSES=0 means: use all. 
NPROCESSES = os.cpu_count() - 2 
########################################################################


ORDER7_FILENAME = "Hurwitz_order7.txt"


def find_type_3B():
    for i in range(1000):
        a0 = Xsp2_Co1('r', 'G_x0')
        o = a0.order()
        if o % 3 == 0:
            a = MM0(a0 ** (o//3))
            if a.chi_G_x0()[0] == 53:
                if Xsp2_Co1(a).subtype[1] in [0,8]:
                    a = MM0(Xsp2_Co1(a))
                    #print(Xsp2_Co1(a), a.order(), 
                    #    a.chi_G_x0(), Xsp2_Co1(a).subtype
                    #)
                    return a
    raise ValueError("No type-3B element found")


g3 = find_type_3B()
#print("g3 =", g3)

z = MM0('x', 0x1000)



def mm_rand(quality = 7):
    l = "pdyx" + "tpl" * quality
    data = [ (t, 'n') for t in l]
    return MM0(data)



def display_mm_rand_data():   
    for i in range(10):
        g = mm_rand()
        g1 = MM(g).reduce()
        print(g1.order(), g1, "\n")

#display_mm_rand_data()




def mm_rand_collect(n = 1000):
    a = np.zeros(120, dtype = np.uint32)
    for i in range(n):
        a[mm_rand().order()] += 1
    return a

#print(mm_rand_collect().reshape((-1,10)))


V3 = MMV(3)
v = V3("R")


def find_order(v, order, z, g):
    v1 = v * (g * z) ** order
    if v1.hash() != v.hash():
        return False
    return v1 == v

    
     



def find_elem(ntrials, order = 7, quality = 6, verbose = 0):
    data = []
    t0 = time.process_time()
    for i in range(ntrials):
        ge =  mm_rand(quality)
        g = g3**ge
        if not find_order(v, order, z, g):
            continue
        if (z * g).order() == order:
             if verbose:
                 print("g3 =", g3)
                 print("ge =", ge)
                 print("z * g3**ge has order", oorder)
             data.append((str(g3), str(ge), order))
    t = time.process_time() - t0
    if verbose:
        print("Run time per case is %.3f ms" % (t1 * 1000 / ntrials))
    return data


def parallel_find_elem(args):
    ntrials, order, quality = args
    return find_elem(ntrials, order, quality)


def mp_find_elem(
        ntasks, 
        ntrials, 
        processes = 0, 
        order = 7,
        quality = 6):
    t0 = time.time()
    table_args = [(ntrials, order, quality)] * ntasks
    n_cpu = os.cpu_count()
    processes = min(processes, n_cpu) if processes else n_cpu
    processes = min(processes, len(table_args))
    s = "Python reports %s CPUs on this computer, using %d of them"
    print(s % (n_cpu, processes))
    with Pool(processes = processes) as pool:
        f = parallel_find_elem
        results = sum(pool.map(f, table_args, chunksize = 1), [])
    pool.join()
    t = time.time() - t0
    total_trials = ntasks * ntrials
    print("Total run time: %.4f seconds" % t)
    print("per trial:       %.4f seconds" % (t / total_trials))
    print("%d trials, %d elements found" % (total_trials, len(results)))
    try:
        f = open(ORDER7_FILENAME, "at")
    except:
        f = open(ORDER7_FILENAME, "wt")
    now_ = datetime.now().isoformat()
    print("# Time: ", now_, ", duration %.2f seconds" % t, file = f )
    print("# %d cases tested, %d cases found" % 
         (total_trials, len(results)), file= f )
    for g3, ge, o in results:
        print("g3 ='%s'" % g3, file = f)
        print("ge ='%s'" % ge, file = f)
        print("order =", o, file = f)
    f.close()
    return results




if __name__ == "__main__":
    #freeze_support()
    order = 7
    while True:
        res = mp_find_elem(100, 500, NPROCESSES, order, 7)
        for g3, ge, o in res:
            print("g3 =", g3)
            print("ge =", ge)
            print("order =", o)





