"""The module contains class ``BitMatrixMulFix``.

Class ``BitMatrixMulFix`` contains a directive that generates
code to multiply a fixed bit matrix with a variable bit matrix.

The C function ``mat24_perm_to_matrix`` uses this kind of matrix
multiplication to convert a permutation in the Mathieu group
to a bit matrix operating on a Golay code word.
"""

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

from random import randint
import types
import sys
import re
import os
from operator import __or__, __xor__
from functools import reduce

from mmgroup.bitfunctions import bitlen, bitweight, hibit, v2
from mmgroup.generate_c import UserDirective


def show_addition_tree(dict):
    print("Addition tree used in module make_addition_table") 
    for x in dict:
        print( hex(x), list(map(hex, dict[x])) )
    print( "" )
    



def small_addition_tree(data):
    data1 = []
    for x in data:
        if bitweight(x) > 1 and not x in data1:
             data1.append(x)
    if len(data1) == 0:
        return {}
    data1 = sorted(data1, key = lambda x: (-bitweight(x),-x))
    x = data1[0]
    for y in data1[1:]:
        if x & y == y:
            data1[0] = x ^ y
            dict = small_addition_tree(data1[:])
            dict[x] = (y, x ^ y)
            return dict
    w, y1 = 1, None
    for y in data1[1:]:
        w1 = bitweight(x & y)
        if w1 > w:
            w, y1 = w1, y
    if not y1 is None: 
        data1[0] = z =  x & y1
        data1 += [x ^ z, y1 ^ z]
        dict = small_addition_tree(data1[:])
        dict[x] = (z, x ^ z)
        dict[y1] = (z, y1 ^ z)
        return dict
    y =  1 << hibit(x)
    data1[0] ^= y
    dict = small_addition_tree(data1[:])
    dict[x] = (y, x^y)  
    return dict  


def make_addition_tree(data, granularity = 8):
    dict = {}
    if len(data) == 0 or max(map(bitweight, data)) <= 1:
        return dict
    maxbit = bitlen(  reduce(__or__, data, 0) )
    stage = {}
    st = 0
    for i in range(0, maxbit, granularity):
        mask1 = (1 << (i + granularity)) -  (1 << i)
        mask  = (1 << (i + granularity)) -  1
        a = [] 
        for x in data:
            xm1 = x & mask1
            xm  = x & mask 
            xm0 = xm1 ^ xm
            if bitweight(xm0) and bitweight(xm1): 
                dict[xm] = (xm1, xm0)
                stage[xm] = st
            
            if bitweight(xm1) > 1 and not xm1 in a:
                a.append(xm1)
        dict1 =  small_addition_tree(a)
        for d in dict1.keys():
            dict[d] =  dict1[d]
            stage[d] = st+1
        st += 2
    #show_addition_tree(dict) 
    for x in (dict):
        for y in dict[x]:
            if bitweight(y) > 1: assert y in dict.keys(), (hex(x),hex(y))

    return dict, stage

       
    
def calc_fanout(data, dict):
    big = 2*len(dict) + 2
    fanout = {}
    for x in dict.keys():
        fanout[x] = big if x in data else 0
    for x in dict.keys():
        for y in dict[x]:
            if y in dict.keys():
                 fanout[y] += 1
    return fanout
     


class register_set(object):
    def __init__(self):
        self.free = []
        self.n = 0
    def alloc(self):
        if not len(self.free):
            self.free.append(self.n)
            self.n += 1
        return self.free.pop()
    def dealloc(self,n):
        assert 0 <= n < self.n
        self.free.append(n)
    def alloc_index(self,n):
        while self.n <= n:
           self.free.append(self.n)
           self.n += 1
        i = self.free.index(n)
        del self.free[i]
        return n
                   
    


def reassign_registers(data, registers):
    mapping = {}
    free_registers =  register_set()
    for i,x in enumerate(data):
        if x in registers.keys() and not registers[x] in mapping.keys():
             mapping[registers[x]] = free_registers.alloc_index(i)
    for r in registers.values():
        if not r in mapping.keys():
             mapping[r] = free_registers.alloc()
    new_registers = {}
    for x in registers.keys():
        new_registers[x] = mapping[registers[x]]
    return new_registers

       



def make_addition_operations(data, granularity = 4):
    calc,stage =  make_addition_tree(data, granularity)
    fanout =  calc_fanout(data, calc)   
    done = set()
    remain = set(calc.keys())
    big = 2*len(calc) + 2
    registers = {}
    ops = []  
    allocated_regs = register_set()  

    def avail(x):
        return bitweight(x) <=1  or x in done   


    while len(remain):
        fo, y1 = [big*big,big*big]  , None
        x = None
        for y in remain:
             doable = avail(calc[y][0]) and avail(calc[y][1])
             fo_new = [stage[y],fanout[y]] 
             for predec in  calc[y]:
                 if predec in calc.keys() and fanout[predec] == 1:
                     fo_new[1] -= 1
             if doable and fo_new < fo:
                 fo, y1 = fo_new, y                

        for predec in calc[y1]:
            if predec in calc.keys():
                 fanout[predec] -= 1
                 if fanout[predec] == 0:
                     allocated_regs.dealloc(registers[predec])
        registers[y1] =  allocated_regs.alloc()
        ops.append((y1,)+calc[y1])
        remain.remove(y1)
        done.add(y1)
    registers = reassign_registers(data, registers)
    return  ops, registers


def do_test_addition_operations(data, ops, registers):
    def getreg(y):
        if bitweight(y) == 1: return y
        else: return results[registers[y]]
    data = list(data)

    ## test matrix calculation

    results = [None] * (max(registers.values()) + 1)
    for x in ops:
         results[registers[x[0]]] = getreg(x[1]) ^   getreg(x[2])
    for i, x in enumerate(data): 
        if bitweight(x) > 1:
            assert x ==  results[registers[x]], map(hex,[x,results[registers[x]]]) 
        if data.index(x)  == i:
            assert registers[x] == i





class BitMatrixMulFix(object):
    """Models multiplication of a fixed by a variable bit matrix.


    The left factor is a fixed bit matrix. Bit matrices are encoded
    as arrays of integers as in module bitfunctions.

    The purpose of this class is to generate C code for perfroming
    this multiplication.
    """

    def __init__(self, granularity = 4):
        """Here 'fixed_matrix' is the fixed left bit matrix factor. 

        'granularity' means something like the number of columns of
        the fixed matrix that are processed in one stage. The default
        appears to be appropriate for m times n matrices,
        12 <= m, n <= 24. 
        """
        self.granularity = granularity

    def set_fixed_matrix(self, fixed_matrix):
        """Here 'fixed_matrix' is the fixed left bit matrix factor. 

        Here we use some magic (vulgo: poorly documented) method
        to obtain a straight-line program calculating the
        non-trivial lines of the matrix product via xor operations.

        If you just want to understand the result (not the magic),
        it is best to check member function selftest.

        """
        self.data = list(fixed_matrix)
        self.nrows =  len(self.data)  
        self.ncols = bitlen(reduce(__or__, self.data, 0))
        self.ops, self.registers = make_addition_operations(
             fixed_matrix, self.granularity)
        self.n_registers = max(self.registers.values()) + 1
        self.n_registers = max(self.n_registers, self.nrows)
        self.copylist = self._make_copylist()

    def _make_copylist(self):
        copylist = []
        for i, x in enumerate(self.data):
            if not x in self.registers.keys() or self.registers[x] != i:
                copylist.append( (i,x) )
        return copylist

              
    def selftest(self):
        """test matrix multiplication fixed_matrix *  unit_matrix


        """
        def getreg(y):
            if bitweight(y) <= 1: return y
            else: return results[self.registers[y]]

        ## test addition operation
        do_test_addition_operations(self.data, self.ops, self.registers)

        ## test matrix calculation

        results = [None] * self.n_registers
        for x in self.ops:
            results[self.registers[x[0]]] = getreg(x[1]) ^ getreg(x[2])
        for i, x in self.copylist:
            results[i] = getreg(x)
        for i, x in enumerate(self.data):
            assert results[i] == x



    def set_names(self, type_name, input_name, output_name):
        self.type_name = type_name
        self.input_name = input_name
        self.output_name = output_name
        self.tmp_name = "_r"
        i = 0
        while self.tmp_name in [self.input_name, self.output_name]:
            self.tmp_name = "_r"+str(i)
            i += 1

    def array_entry(self, array_name, index):
        return array_name + "[" + str(index) + "]"

    def var_name(self, x):
        if x & (x-1):
            r = self.registers[x]
            if r >= self.nrows:
                return self.array_entry(self.tmp_name, r-self.nrows)
            else: 
                return self.array_entry(self.output_name, r)
        elif x:
            r = v2(x)
            assert x == 1 << r
            return "(%s)(%s)" % (self.type_name,
                                      self.array_entry(self.input_name, r))
        else:
            return "0"

    def make_c_calculation(self):
        s = ""
        for ops in self.ops:
            s +=  " %s = %s ^ %s;\n" % tuple(map(self.var_name, ops))
        for dest, src in self. copylist:
            s += " %s = %s;\n" % (
               self.array_entry(self.output_name,dest), self.var_name(src))
        return s;

    def make_c_comment(self):
        s = """GF(2) bit matrix multiplication as a sequence of XOR operations:
  %s  = _MFIX_ * %s  ,
with bit matrices represented as arrays of unsigned integers,
each integer representing a row of the matrix. 
Here _MFIX_ is the fixed bit matrix  {\n""" % (
   self.output_name, self.input_name)

        w = max(1, (self.ncols + 3) >> 2)
        l = (60) // (w+3)
        for i, d in enumerate(self.data):
            if i%l == 0:  s += " "
            s += "0x%0*x" % (w,d)
            if i ==  len(self.data)-1: 
                s += "\n}"
            elif i%l == l-1: 
                s += ",\n"
            else: 
                 s += ","
        lines = [  "// " + y for y in s.split("\n")] 
        return "\n".join(lines)   + "\n"

    def make_c_program(self):
        s = "\n" + self. make_c_comment() +  "{" + "\n"
        if self.n_registers > self.nrows:
            s +=  " %s %s[%d];\n" % (self.type_name,self.tmp_name,
                 self.n_registers - self.nrows)
        s += self.make_c_calculation()
        s += " // %d operations generated\n" % len(self.ops)
        s +=  "}" + "\n"
        return "\n" + s + "\n"
        
                        
 
    def generate_c_bitmatmul(self, table, type, input_var, output_var): 
        self.set_fixed_matrix(table)
        self.set_names(type, input_var, output_var)
        self.selftest()
        return self.make_c_program()
 

    def tables(self):
        """This class generates no tables"""
        return {}

    def directives(self):
        """Code generation directive has name "BITMATMUL". """
        return {
           "BITMATMUL" : UserDirective(self.generate_c_bitmatmul, "psss"),
        }
    
   