from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import random
from random import randint
import types
import sys
import re
import os
from operator import __or__, __xor__
import subprocess

import pytest
 
from mmgroup.tests.tempdir import make_temp_dir, kill_temp_dir
from mmgroup.tests.compile import compile_testprogramm
from mmgroup.dev.mat24.make_mul_transp import BitMatrixMulTransp
#from mmgroup.dev.mat24.make_mul_transp import generate_c_mul_transp
  

#########################################################################
######### Tests #########################################################
#########################################################################



@pytest.mark.slow
@pytest.mark.mat24
@pytest.mark.compiler
def test_make_mul_transp():
    """Generate, compile and test a C program"""
    tmp_dir =  make_temp_dir()
    FILE_NAME = "tmp_test_mul_transp"
    NAME = os.path.join(tmp_dir, FILE_NAME)
    C_NAME, EXE_NAME = NAME + ".c", NAME + ".exe"
  
    def make_testprogram(lv, lm, n=1):
        generate = BitMatrixMulTransp().generate_c_mul_transp
        f = open(C_NAME, "wt")
        print(r"""#include <stdio.h>
#include <stdint.h>
int main(int argc, char **argv)
{{
    uint_fast32_t i, v;
    static uint_fast32_t a[1000];
    sscanf(argv[1], "%ld", &v);
    for (i=0; i < argc-2; ++i) sscanf(argv[i+2], "%ld", a+i);
    {0}
    printf("%ld\n",v);
    return 0;
}}
""".format(generate("v", "a", lv, lm, n)), file = f)
        f.close()
        #subprocess.check_output(["gcc", C_NAME, "-o", EXE_NAME])
        compile_testprogramm([C_NAME], EXE_NAME)
        checker = BitMatrixMulTransp()
        checker.set_matrix(lv, lm, n)
        return checker

    def run_testprogram(v, m):
        data = list(map(str, [v]+m))  
        res =  subprocess.check_output([EXE_NAME] + data)
        return int(res)  

    def test_testprogram(checker , v, m):
        res =  run_testprogram(v, m)
        ref = checker.compute(v,m+[0]*64)
        assert ref == ref, (hex, ref)

    def del_testprogram():
        for f in [C_NAME, EXE_NAME]:
            os.remove(f)

    test_data_dict = {
        (4,4,1,0): [ (3, [2,3, 1, 1]) , ]
    }

    def test_data(lv, lm, n, n_tests=10):
        """yield  n_tests test cases (v, m) for given lv, lm, m0"""
        try:
            d =  test_data_dict[(lv, lm, n, m0)]
            for v, m in d: yield v, m
        except:
            pass
        for i in range(n_tests):
            v = random.randint(0, 2**(n*lv)-1)
            m = [] 
            for i in range(lm):
                 m.append(random.randint(0, 2**(n*lv)-1))
            yield v, m

    def test_programs():
        """yield cases (lv, lm, n) for a making a test program"""
        yield 4, 4, 1
        yield 16, 12, 2

    

       
    for lv, lm, n in test_programs():
        print("Make C program for case lv=%d, lm=%d, n=%d" % (
                 lv,lm,n))
        checker = make_testprogram(lv, lm, n)  
        print ("C program has been made, starting tests..")
        for v, m in test_data(lv, lm, n):
              test_testprogram(checker, v, m)
        del_testprogram()
        print ("passed")

    kill_temp_dir(tmp_dir)


