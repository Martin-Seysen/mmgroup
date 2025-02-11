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



#######################################################################
# Auxiliary functions for function ``make_addition_tree``
#######################################################################

def small_addition_tree(data):
    """Auxiliary function for function ``make_addition_tree``

    The function returns a dictionary with entries of shape 

    ``r:  (x, y)``

    for integers ``r, x, y`` such that ``r = x ^ y``, and  each of 
    the entries ``x, y`` is either of bit weight 1 or it is also
    contained as a key in the dictionary.

    The dictionary contains all entries of the list ``data`` of
    bit weight > 1 as keys.

    This function is highly recursive and should be used only if
    the bit weight of the union of the entries of data
    is not too large (say, at most 8).
    """
    # Collect all unprocessed entries of array ``data`` in ``data1``
    # with bit weight at least 2
    data1 = []   # list of entries yet to be processed
    for x in data:
        if bitweight(x) > 1 and not x in data1:
             data1.append(int(x))
    # Done if the list ``data1`` is empty.
    if len(data1) == 0:
        return {}
    # Sort list ``data1`` by descending bit weight 
    data1 = sorted(data1, key = lambda x: (-bitweight(x),-x))
    # Process first entry ``x`` of data1
    x = data1[0]
    # If list ``data1`` contains a word ``y`` with ``y = x & y``
    # then append an entry ``x: (y, x ^ y)`` to the output dict
    # and replace entry ``x`` in ``data1`` by ``x ^ y``.
    # Then recursively enter the result of
    # ``small_addition_tree(data1)`` to output dictionary ``dict``
    # and return.
    for y in data1[1:]:
        if x & y == y:
            data1[0] = x ^ y
            dict = small_addition_tree(data1[:])
            dict[x] = (y, x ^ y)
            return dict
    # If the previous attempt has failed then find an entry ``y1``
    # (different from ``x```) in the list such that the bit 
    # weight of ``x & y1`` is maximal.
    w, y1 = 1, None
    for y in data1[1:]:
        w1 = bitweight(x & y)
        if w1 > w:
            w, y1 = w1, y
    # If a suitable ``y1`` of bit weight at least 2 has been found
    # then delete entry ``x`` from the list ``data1`` add the entry 
    # ``z`` with ``z = x & y1`` and also the entries ``x ^ z``, 
    # ``y1 ^ z``  to the list. The let the output dict consist of the 
    # entries ``x: (z, x ^ z),  y1:(z, y1 ^ z)`` and the recursively
    # computed output ``small_addition_tree(data1)`` and return
    if not y1 is None: 
        data1[0] = z =  x & y1
        data1 += [x ^ z, y1 ^ z]
        dict = small_addition_tree(data1[:])
        dict[x] = (z, x ^ z)
        dict[y1] = (z, y1 ^ z)
        return dict
    # If no suitable ``y1`` of bit weight at least 2 has been found
    # then let ``y`` be the power of two with the highest bit of ``x``
    # being set. Then append an entry ``x: (y, x ^ y)`` to the output 
    # and and replace entry ``x`` in ``data1`` by ``x ^ y``.
    # Then recursively enter the result of
    # ``small_addition_tree(data1)`` to output dictionary ``dict``
    # and return.
    y =  1 << hibit(x)
    data1[0] ^= y
    dict = small_addition_tree(data1[:])
    dict[x] = (y, x^y)  
    return dict  


#######################################################################
# Function ``make_addition_tree``
#######################################################################


def make_addition_tree(data, singleton = False, granularity = 8):
    """Create a direct addition tree for a set of integers

    The function returns a pair ``(dict, stage)``. Here 
    dictionary ``dict`` has entries of shape ``r:  (x, y)``
    for integers ``r, x, y`` such that ``r = x ^ y``, and  each of 
    the entries ``x, y`` is either of bit weight 1 or it is also
    contained as a key in the dictionary.

    The dictionary contains all entries of the list ``data`` of
    bit weight > 1 as keys.

    This function calls ``small_addition_tree`` (possibly several 
    times) with a list of data, such that the bit weight of 
    the union of these data is at most equal to ``granularity``.

    if parameter ``singleton`` is True then all bit positions 
    set in any entry of ``date`` are entered into the dictionary 
    as key of bit weight 1, with is value being an empty tuple.

    Dictionary ``shape`` has the same set of keys as dictionary
    ``dict``. For each key the value in that dictionary indicates 
    a **stage** at which this entry is computed.
    """
    dict = {}
    stage = {}
    if len(data) == 0 or max(map(bitweight, data)) <= 1:
        return dict, stage
    all_bits =  reduce(__or__, data, 0)
    maxbit = bitlen(all_bits)
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
                stage[xm] = st+1
            
            if bitweight(xm1) > 1 and not xm1 in a:
                a.append(xm1)
        dict1 =  small_addition_tree(a)
        for d in dict1.keys():
            dict[d] =  dict1[d]
            stage[d] = st
        st += 2
    # Check addition tree in ``dict``
    #show_addition_tree(dict) 
    for x in (dict):
        for y in dict[x]:
            if bitweight(y) > 1: assert y in dict.keys(), (hex(x),hex(y))
    # Add singeltons to ``dict`` if requested
    if singleton:
        while all_bits:
            entry = all_bits & -all_bits
            dict[entry] = ()
            all_bits &= ~entry

    return dict, stage

#######################################################################
# Auxiliary functions for function ``make_addition_operations``
#######################################################################
      
    
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
    def dealloc(self, n):
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

       



def make_add_operations(data, singleton = False, granularity = 4):
    calc,stage =  make_addition_tree(data, singleton, granularity)
    fanout =  calc_fanout(data, calc)   
    done = set()
    big = 2*len(calc) + 2
    registers = {}
    ops = []  
    allocated_regs = register_set() 

    for x, d in calc.items():
        if bitweight(x) <= 1:
            registers[x] = allocated_regs.alloc()
            done.add(x) 
            ops.append((x, x))
        for y in d:
            if bitweight(y) <= 1:
                done.add(y)
    remain = set(calc.keys()) - done


    while len(remain):
        fo, y1 = [3, big*big]  , None
        x = None
        for y in remain:
             doable = calc[y][0] in done and calc[y][1] in done   
             fo_new = [2, stage[y]] 
             for predec in  calc[y]:
                 if predec in calc.keys() and fanout[predec] == 1:
                     fo_new[0] -= 1
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


def copy_operation(data, ops, registers):
    copies = []
    zeros = []
    for i, x in enumerate(data):
        if not x in registers.keys() or registers[x] != i:
           if x:
                copies.append( (i, x, x) )
           else:
                zeros.append(i)
    for i in zeros:
        copies.append( (i, 0, 0) )
    return copies


#######################################################################
# Function ``make_addition_operations``
#######################################################################


def make_addition_operations(m, singleton = False, granularity = 4):
    """Straight line program multiplying fixed with variable bit matrix

    The function returns a staight-line program for multiplying
    a fixed bit matrix ``m`` with a variable bit matrix ``a``. Here
    a bit matrix is encoded as an array (or a list) of integers as in 
    module ``mmgroup.bitfunctions``. Matrix ``m`` is the left factor 
    of the matrix product. Multiplication with matrix ``m`` from
    the left corresponds to a sequence of row operations applied to
    matrix ``a``.

    The straight line program for computing ``m * a`` is returned
    as a array of tuples ``(register, value, op1 [, op2])``, where
    each tuple encodes a binary or a copy operation. If entry 
    ``op2`` is not present then we have a copy operation. In a 
    nutshell this tuple encodes the operation:

    Compute ``(value * a) = (op1 * a) + (op2 * a)`` and store the 
    result in register ``register``.

    ``register`` is the number of register where to store the result 
    of the operation. Register numbers start with 0. Eventually,
    ``register[i]`` will contain the ``i``-th row of the output
    matrix. That number of rows of the result ``m * a`` is the 
    maximum of the bit lengths of the entries of input entry ``m``; 
    so the last row of the output matrix is the highest column 
    of ``m`` with nonzero  entries. More registers may be required 
    for storing intermediate results.
   
    Operands and results of an operation are bit vectors corresponding
    to rows of bit matrices. Here the value ``v`` of an entry 
    ``value``, ``op1``, or ``op2`` is an integer encoding a bit vector 
    ``v``, such  that this entry is equal to the bit vector ``v * a``. 
    If ``op2`` is present then the operation encoded by this tuple is 
    the addition of these two operands (as bit vectors); and we 
    have ``value = op1 ^ op2``. Otherwise we have a copy operation
    an then we have  ``value = op1``.

    In case ``v = 0`` the corresponding bit vector is zero; and if
    ``v`` has bit weight 1 then the corresponding bit vector is a row
    of matrix ``a``.  A bit vector corresponding to any other value 
    of ``v`` must be (and will be) computed in a register.

    So the user may generate e.g. C code for computing ``m * a``
    from ``a`` by using the output of this function. Here the user
    should maintain a dictionary that maps a value ``v``
    (corresponding to a bit vector ``v * a``) to the number of
    the register containing that bit vector. This dictionary
    is required for locating non-trivial operands in the 
    register set.

    If parameter ``singleton`` is True then all rows of the input
    matrix ``a`` are loaded into the registers before any binary
    operations will be done. 

    Parameter ``granularity`` controls the construction of the
    addition tree as in function ``make_addition_tree``.
    """
    ops, regs = make_add_operations(m, singleton)
    copies = copy_operation(m, ops, regs)
    ops_out = []
    for op in ops:
        d = (regs[op[0]],) + op
        ops_out.append(d)
    return ops_out + copies



#######################################################################
# Testing function ``make_addition_operations``
#######################################################################



def check_addition_operations(m, ops, singleton = False):
    """Check result of function ``make_addition_operations``

    Here parameters ``m`` and ``singletons`` are inputs to
    function ``make_addition_operations`` and parameter ``ops``
    is the corresponding return value of that function.

    Then this function ``check_addition_operations`` performs a
    straight-line program computing ``m * a`` for the unit  matrix 
    ``a`` based on parameter ``ops``. It raises ValueError if the
    result of that straight-line program differs from ``m``.

    If parameter ``singleton`` is True than we also check the
    condition for this case as stated in the documentation of
    function ``make_addition_operations``.    
    """
    num_regs = max([op[0] for op in ops]) + 1
    results = {}
    registers = {}
    
    def check_operand(value):
        in_registers = value in registers
        if  bitweight(value) > 1:
            assert in_registers
        if in_registers:
            assert results[registers[value]] == value
        return in_registers

    for op in ops:
        reg, value, operands = op[0], op[1], op[2:]
        assert len(operands) > 0
        in_registers = min([check_operand(value) for value in operands])
        if len(operands) > 1 and singleton:
            assert in_registers
        assert value == reduce(__xor__, operands)
        results[reg] = value
        registers[value] = reg

    for i, x in enumerate(m):
        assert results[i] == x
    return


#######################################################################
# Display result of function ``make_addition_operations``
#######################################################################


def display_addition_operations(ops):
    """Display result of function ``make_addition_operations``

    Let ``ops`` be a result of a call to function 
    ``make_addition_operations``. Then we display the straight-
    line program corresponding to that result in a readable form.
    """
    print("Straight-line program generated by make_addition_operations") 

    registers = {}
    
    def display_operand(value):
        if value in registers:
            return 'o[%d]' % registers[value]
        if value == 0:
            return "0"
        if  bitweight(value) == 1:
            return "i[%d]" %  v2(value)
        raise ValueError("Illegal operand in straight line program")

    for op in ops:
        reg, value, operands = op[0], op[1], op[2:]
        op_str = " ^ ".join([display_operand(s) for s in operands])
        assign_str = "%-5s = %s" % ("o[%d]"  % reg, op_str)
        print("%-30s // 0x%08x" % (assign_str, value))
        registers[value] = reg


#show_addition_tree = display_addition_operation

#######################################################################
# Class ``BitMatrixMulFix`` wrapping ``make_addition_operations``
#######################################################################


class BitMatrixMulFix(object):
    """Models multiplication of a fixed by a variable bit matrix.


    The left factor is a fixed bit matrix. Bit matrices are encoded
    as arrays of integers as in module bitfunctions.

    The purpose of this class is to generate C code for performing
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

    def set_fixed_matrix(self, fixed_matrix, singleton = False):
        """Here 'fixed_matrix' is the fixed left bit matrix factor. 

        Here we use some magic (vulgo: poorly documented) method
        to obtain a straight-line program calculating the
        non-trivial lines of the matrix product via xor operations.

        If you just want to understand the result (not the magic),
        it is best to check member function selftest.

        """
        self.data = list(fixed_matrix)
        self.singleton = singleton
        self.nrows =  len(self.data)  
        self.ncols = bitlen(reduce(__or__, self.data, 0))
        self.ops = make_addition_operations(
            fixed_matrix, singleton, self.granularity
        )
        self.n_registers = max([op[0] for op in self.ops]) + 1
        self.selftest()



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

    def register_name(self, reg):
        if reg >= self.nrows:
            return self.array_entry(self.tmp_name, reg - self.nrows)
        else: 
            return self.array_entry(self.output_name, reg)

    def operand_name(self, x, registers):
        if x in registers:
            return self.register_name(registers[x])
        if x == 0:
            return "0"
        r = v2(x)
        assert x == 1 << r
        return "(%s)(%s)" % (self.type_name,
                                self.array_entry(self.input_name, r))
 
    def make_c_calculation(self):
        s = ""
        registers = {}
        for op in self.ops:
            reg, value, operations = op[0], op[1], op[2:]
            res_name = self.register_name(reg)
            operand_names = [
                self.operand_name(x, registers) for x in operations]
            registers[value] = reg
            operation = " ^ ".join(operand_names)
            s +=  " %s = %s;\n" % (res_name, operation)
        return s;

    def ops(self):
        return self.ops

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
        return self.make_c_program()

    def selftest(self):
        check_addition_operations(self.data, self.ops, self.singleton)
 

    def tables(self):
        """This class generates no tables"""
        return {}

    def directives(self):
        """Code generation directive has name "BITMATMUL". """
        return {
           "BITMATMUL" : UserDirective(self.generate_c_bitmatmul, "psss"),
        }
    
   