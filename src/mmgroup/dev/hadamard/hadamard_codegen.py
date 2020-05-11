"""Generate C code for Hadamard-like matrixes.

In [Seys19] we define two generators tau and xi of the monster
group and their representation matrices. Most of the non-monomial
parts of these matrices can be implemented as Hadamard matrices, 
which must be scaled by negative power of two. In order to obtain
tau and xi, the Hadamard matrices must also be multiplied by 
certain monomial matrices from the left and from the right.

Here a Hadamard matrix is a 2**n times 2**n matrix H with indices
running from 0 to 2**n - 1. Matrix H has entries 
H[i,j] = parity(i & j), with the '&' operator as in C and
parity(x) the bit parity of the nonnegative integer x, when x is 
given in binary notation.

This module contains a class HadamardMatrixCode for generating
C code that multiples a vector v with a Hadamard matrix H from 
the right. The entries of the vector v are integers modulo p for 
a fixed (not necessarily prime) number p = 2**k-1, 2 <= k <= 8.
 
Vector v is coded as an array of integers of type uint_mmv_t,
where uint_mmv_t is uint32_t or uint64_t, depending on the
settings in module config.py. Here an single integer of type 
uint_mmv_t may store several components of vector v as described
in module mm/mm_aux.py. Class HadamardMatrixCode is a subclass
of class mm/mm_basics.MM_Op, which contains basic operations for
dealing with such a vector.

Our hope is that an optimzing C compiler completely performs the 
relevant operations in registers, if sufficiently many of them 
are present. Note that fewer operations are necessary if v consists 
of 64-bit integers. On the other hand, modelling 64-bit integers on 
a 32-bit architecture causes unnecessary overhead. So we support 
both, 32-bit and 64-bit integers.  



References
----------
see file references.txt
"""



#23456789012345678901234567890123456789012345678901234567890123456789012345

from __future__ import absolute_import, division, print_function
#from __future__ import  unicode_literals



import sys
import os
from numbers import Integral

from mmgroup.bitfunctions import bitparity, bitweight
from mmgroup.dev.mm_op.mm_op import MM_Op, INT_BITS, c_hex



def as_c_expr(x):
    """Map integers to hex string, don't change anything else

    This function is used for entering objects into automatically
    generated C programs in the following class C_Expr, which 
    models expressions in the C language.  

    Most of these objects have a 'str' or '__repr__' method that
    formats the object as required for the C language.

    This would also work for integer constants, when writing 
    them as decimal numbers. In the following class 
    HadamardMatrixCode almost all integer constants are bit
    masks, so we better write them as hexadecimal constants
    of sufficient bit length. Function c_hex in module mm_basics 
    does this job.
    """
    if isinstance(x, Integral):
        return c_hex(x)
    return x



class C_Expr:
    """Models an arbitrary C expression

    Applying a standard operator such as '+', '-', '&', '|', '^', 
    '>>', '<<',  to an instance of this class returns the C code 
    for performing the corresponding operation as a string. 

    Atomic C expressions, e.g. variables, must be given as instan-
    ces of subclasses of this class; see e.g. class C_UintVar.
    One operand of an operation with an instance of this class
    may be an integer. Then the corresponding integer constant
    (written hexadecimally) is used as that operand. Note that
    shift factors in the '<<' and '>>' operators must be integers.

    Here only the C operators required for our purposes are 
    implemented. We do not make any attempt to save brackets.
    """
    def __init__(self, name):
         self.name = name

    def __repr__(self):
        return self.name

    str = __repr__

    def __lshift__(self, other):
         assert isinstance(other,int), "2nd shift operand must be integer"
         if other == 0:
              return self
         elif other > 0:
              return C_Expr("(%s << %d)" % (self, other) ) 
         elif other < 0:
              return C_Expr("(%s >> %d)" % (self, -other) ) 

    def __rshift__(self, other):
         return self.__lshift__(-other)

    def _binop(self, other, op):
         oper1 = str(as_c_expr(self))
         oper2 = str(as_c_expr(other)).strip()
         if len(oper1) + len(oper2) > 40:
             op_ = "\n    " + op + " "
         else:
             op_ = " " + op + " "
         return  C_Expr("(" + oper1 + op_ + oper2 + ")") 

    def __add__(self, other):
        return self._binop(other, '+')

    __radd__ = __add__

    def __sub__(self, other):
        return self._binop(other, '-')

    def __and__(self, other):
        return self._binop(other, '&')

    __rand__ = __and__

    def __or__(self, other):
        return self._binop(other, '|')

    __ror__ = __or__

    def __xor__(self, other):
        return self._binop(other, '^')

    __rxor__ = __xor__

    def __invert__(self):
        return  C_Expr("(~%s)" % self)

    def __neg__(self):
        return  C_Expr("(-%s)" % self)

    def assign(self, value):
        """Return string encoding an assignment statement: self = value;
     
        'value' must be another instance of this class or a string  
        containing a valid C expression or an integer.
        """ 
        return  "%s = %s;\n" % (self, as_c_expr(value))
 
    def assign_xor(self, value):
        """Return string encoding an assignment statement: self ^= value;
     
        'value' must be another instance of this class or a string  
        containing a valid C expression or an integer.
        """ 
        return  "%s ^= %s;\n" % (self, as_c_expr(value))

    def assign_and(self, value):
        """Returnstring encoding an assignment statement: self &= value;
     
        'value' must be another instance of this class or a string  
        containing a valid C expression or an integer.
        """ 
        return  "%s &= %s;\n" % (self, as_c_expr(value))

    def assign_or(self, value):
        """Return string encoding an assignment statement: self |= value;
     
        'value' must be another instance of this class or a string  
        containing a valid C expression or an integer.
        """ 
        return  "%s |= %s;\n" % (self, as_c_expr(value))


class C_UintVar(C_Expr):
    """Models a special integer variable of a given type

    Instances of this class model integer variables which are designed 
    for use with C operations as specified in class C_Expr.

    Here a variable has a name and a type, which must be a valid
    integer type in the C code to be generated.

    Such variables are also designed to be a member of a list of 
    variables. Each variable as an optional attribute 'index'
    specifying the index of the varible in that list.
    """
    def __init__(self, c_type, name, index = None, var_number = None):
         """Create a C variable with a given type and 'name'.

         The type of the variable is the string given in the
         argument 'c_type'. 'index' is an optional argument
         specifiying the index of the variable in a list.
         This list is kept by the code generating  process 
         only; and it is irrelevant for the C compiler.

         Usually variables are created as instances of this
         class when needed. Then their names contain
         consecutive numbers, e.g. 'r0', 'r1', 'r2', ...    
         The corresponding number of a variable may be 
         passed with argument 'var_number'.
         """
         self.name = name
         self.var_number = var_number
         self.index = index
         self.c_type = c_type 


class C_UintVarPool:
    """Model a (usually short) list of integer variables.

    Assume that we want to modify a short vector v of integers 
    stored in memory. Therefore we load all entries of vector v
    into variables of the appropriate type. Then we perform the
    requested operations on these variables. Finally, we store
    all variables back to memory. Thus an optimizing compiler
    may save run time by keeping some or all variables in 
    registers. 

    So we may create a pool of variables of a given type with
    a given name. Here name = "r" means that the variables of
    the pool with have names "r0", "r1", "r2",... .

    All variables in the pool are instances of class C_Expr,
    so we may perfor C operations on them as decribed in 
    class C_Expr.

    The pool has an 'official'  size equal to len(self),
    which is set in the constructor, and which can be modified
    by method self.extend_pool(new_size). The official
    variables are ordered, so self[i] is the  i-th official
    variable in the pool. This mean that the official variables 
    may represent a vector of integers. The order of the
    variables can be modified with method self.xch(); so
    exchanging components of the vector is for free.

    Apart from the official variables there are temporary
    variables of the same type in the pool. self.temp(i) 
    returns the i-th temporary variable. We often need 
    temporary variables for the operations on the entries of 
    the vector. Official variables can also be exchanged
    with temporary variables.

    An instance of this class keeps track of all variables 
    created. Method self.declare() can be used to declare all 
    these variables in a C program.
    """
    LINE_LEN = 50

    def __init__(self, c_type = "uint_mmv_t", name = "r%d", length = 0):
        """Create a pool of variables of type 'ctype'.

        The pool will contain 'length' official variables. These
        variables are instances of class C_UintVar and they are
        obtained as self[0], ..., self[length - 1].

        The names of the variables are name % 0, name % 1, name % 2, ..
        where name is given by parameter 'name'. So for obtaining
        r0, r1, r2, ... parameter 'name' should be "r%d"
        """
        self.a_size = 0      # Number of official variables
        self.vars = []       # List of all veriables, Here 
                             # official variable are listed
                             # in their current order, followed
                             # by the temporary variables.
        self.name = name     # Name string used for variables
        self.c_type = c_type # C type of the variables
        self.resize(length)  # Create 'length' official variables
         

    def extend_pool(self, new_size):
        """Internal method

        The method adds variables to the pool so that there are
        at least 'new_size' variables in the pool.
        """
        while len(self.vars) < new_size:
            l = len(self.vars)
            name = self.name % l
            self.vars.append(C_UintVar(self.c_type, name, l, l))

    def __getitem__(self, i):
        """Return the i-th official variable"""
        assert 0 <= i < self.a_size, (i, self.a_size)
        return self.vars[i]

    def resize(self, new_size):
        """Change the size of the pool to 'new_size'

        Afterwards, the 'offical' part of the pool has has the
        length given by 'new_size'. All temporary variables in
        the pool are discarded.
        """
        self.extend_pool(new_size)
        self.a_size = new_size

    def temp(self, i = 0):
        """Return the i-th temporary variable"""
        self.extend_pool(self.a_size + i + 1)
        return self.vars[self.a_size + i]

    def xch(self, i, j):
        """Exchange self[i] with self[j].

        i and j may also be variables in the pool (including
        temporary variables) insted of their indices. Thus

           v1 = self[1]; t = self.temp(0); self.xch(v1, t)

        exchanges an official variable with a temporary variable.

        This method performs just book keeping, so that we need 
        not generate any code here.   
        """
        if isinstance(i, C_UintVar):
            var_i, i = i, i.index 
        else:
            var_i = self.vars[i]
        if isinstance(j, C_UintVar):
            var_j, j = j, j.index 
        else:
            var_j = self.vars[j]
        assert self.vars[i].index == i
        assert self.vars[j].index == j
        if i == j:
            return
        assert 0 <= i < len(self.vars)
        assert 0 <= j < len(self.vars)
        self.vars[i], self.vars[j] = self.vars[j], self.vars[i]
        self.vars[i].index = i
        self.vars[j].index = j

    def declare(self):
        """Declare all variables in the pool in C

        The function returns a string containing that declaration.
        """
        line_len = self.LINE_LEN - len(self.c_type) - 1
        nvars =  line_len // (len(self.name) + 4)
        pos = 0
        s = ""
        while pos < len(self.vars):
           end_pos = min(pos + nvars, len(self.vars))
           vars = [self.name %i for i in range(pos, end_pos)]
           s += self.c_type + " " + ", ".join(vars) + ";\n"
           pos = end_pos
        return s
         
    def __iter__(self):
        """Iterate over official variables in their current order"""
        for var in self.vars[:self.a_size]:
            yield var 

    def __len__(self):
        """Return the number of the official variables"""
        return self.a_size


    def var_number_list(self):
        """Return list of numbers of the official variables

        Given a 'name' for the variables of the pool, the
        official variables have names <name>0, <name>1, <name>2 ..
        These names may be permuted by method self.xch().

        This method returns the list of numbers in the names
        of the official variables in the current order. 
        """
        return [self.vars[i].var_number for i in range(self.a_size)]


class HadamardMatrixCode(MM_Op):
    """Generate C code for multiplying vector by a Hadamard matrix 
 
    We generate C code for multiplying a vector v of integers
    modulo p = 2**k-1  by a Hadamard matrix H. 

    Vector v is stored in memory as described in class
    mm/mm_aux.py. It is loaded into a set of C variables
    which is constructed as an instance of class 
    C_UintVarPool. C code for multiplication of v by H is 
    generated with the operations given in class C_Expr,
    operating on these variables. Finally the result is
    written back from the variables to memory.

    Matehmatically, a 2**n times 2**n times Hadamard matrix
    is constucted as product of n matrices of shape

       I2  (x) ... (x) I2 (x) H2 (x) (I2) (x) ... (x) (I2) , 
    
    where (x) denotes the Kronecker product, I2 is the 2 times 2 
    unit matrix, and H2 is the 2 times 2 Hadamard matrix. H2 
    occurs exactly once in each of the k factors of the matrix
    product. In different factors, H2 occures at different 
    positions. Note that all these k matrix factors commute. 

    The m-th matrix factor of H is implemented as a butterfly 
    operation on vector v in the usual way:

         v[i], v[i+j] = v[i] + v[i+j], v[i] - v[i+j],

    with j = 2**m, for all integers i where bit m of i is zero. 
    """
    def __init__(self, p, log_vlength, verbose = 0):
        """Create an instance of class c_matrix_code

        Operation on vectors are done modulo p as described
        in class mm/mm_aux.py. p must be of shape
        p = 2**n-1, 2 <= n <= 8.

        'log_vlength' is the binary logarithm of the size of
        the vector v and of the Hadamard matrix H. We test
        the cases 'log_vlength' = 2, 4, 6 corresponding to
        vectors of size 4, 16, 64 only.
        """
        super(HadamardMatrixCode, self).__init__(p)
        self.LOG_VLEN = log_vlength   
        self.NO_CARRY = self.FIELD_BITS == self.P_BITS
        self.verbose = verbose
        self.make_h_masks()
        self.reset_vars()
        if verbose:
           print("HadamardMatrixCode mod", p,  ", size =", self.vlen)


    def reset_vars(self, log_vlen = None):
        """Reset variables 

        This method resets all all information required for the 
        the code generation process, including the pool self.vars 
        of variables, which is an intance of class C_UintVarPool.

        If log_vlen is specified, self.LOG_VLEN is changed to log_vlen.        
        """
        if log_vlen: self.LOG_VLEN = log_vlen
        if self.verbose: print("reset vars")
        self.vlen = max(1, (1 << self.LOG_VLEN) >> self.LOG_INT_FIELDS)
        self.vars = C_UintVarPool("uint_mmv_t", "r%d", self.vlen)
        self.expanded = 0
        self.matrix_code = ""     # C implementation of matrix operation
        self.n_code_lines = 0     # just for staticstics
        self.n_operations = 0     # just for staticstics
      

    def make_h_masks(self):
        """Create bit masks required for a butterfly operation

        self.ALL1_MASK is the all-bits-one mask.

        For a power i of two, self.B0_MASK[i] selects the 
        entries v[j] of a vector stored in a single variable 
        with i & j = 0, and self.B1_MASK[i] selects the 
        entries v[j] with i & j = i.

        self.H_MASK_EXTERN is for butterfly op. on different 
        variables var[i], var[j]. It selects all valid fields 
        of a variable. self.H_MASK_CARRY selects the carry
        bits of all valid fields

        self.H_MASK[i] is for butterfly op. on (v[j], v[j+i]),
        with v[i] nd v[i+j] stored in the same variable, for 
        the relevant powers i of two. Here the mask selects 
        the components v[j] of v with i & j == i.

        Note that in case p = 2**2**k-1, (i.e.  p = 3, 15, 255),
        the bit fields with odd index in a variable are invalid
        when a vector is expanded for a Hadamard operation, as
        explained in method expand_hadamard(). self.H_MASK[j]
        takes care of this special situation.
        """
        self.ALL1_MASK = (1 << self.INT_BITS) - 1
        self.B0_MASK = {} 
        self.B1_MASK = {} 
        self.H_MASK = {} 
        mask = self.smask(self.P, -1, self.FIELD_BITS << self.NO_CARRY) 
        self.H_MASK_EXTERN = mask
        mask = self.smask(self.P+1, -1, self.FIELD_BITS << self.NO_CARRY) 
        self.H_MASK_CARRY = mask
        for i0 in range(self.LOG_INT_FIELDS):
            i = 1 << i0
            mask = [x for x in range(self.INT_FIELDS) if x & i == 0]
            self.B0_MASK[i] = self.smask(self.P, mask)
            mask = [x for x in range(self.INT_FIELDS) if x & i == i]
            self.B1_MASK[i] = self.smask(self.P, mask)
            self.H_MASK[i] = self.B1_MASK[i] & self.H_MASK_EXTERN
       
            
    def comment(self, comment):
        """Enter a 'comment' into the C code"""
        for s in comment.split("\n"):
            self.matrix_code += "// " + s + "\n"

    def comment_vector(self):     
        """Comment of represention of vector v in  variables"""
        vl = self.vars.var_number_list()
        vdata = ",".join(map(str,vl)) 
        if len(vl) > 10: vdata = "\n " + vdata 
        name = self.vars.name.replace("%d", "%s") % "(i)"        
        self.comment("Vector is now  %s for i = %s" % (
            name, vdata))

    def comment_statistics(self, reset = True):
        """Write statistics as comment into the C code"""
        self.comment("%d lines of code, %d operations" %
               (self.n_code_lines, self.n_operations))
        if reset:
            self.n_code_lines = self.n_operations = 0;
     
    def reduce_butterfly(self, var, dest = None):
        """Auxilary function for method butterfly_op().
   
        Performs reduction mod p after a butterfly operation
        for all vomponents v[i] stored in variable 'var'.

        The reduced result is stored in variable 'dest'.
        'dest' defaults to 'var'.
        """
        mask = self.H_MASK_EXTERN 
        if dest is None or dest == var:
            dest, t = var, self.vars.temp()
        else:
            t = dest
        s = t.assign(var & self.H_MASK_CARRY)
        s += dest.assign(var - t + (t >> self.P_BITS))            
        self.matrix_code += s
        self.n_code_lines += 2
        self.n_operations += 4

    def internal_butterfly(self, var, j):
        """Auxilary function for method butterfly_op().
   
        Perform butterfly operation inside the single C variable 
        'var'.  For all components v[i] of v stored in var we put

            v[i], v[i+j] = v[i] + v[j], v[i] - v[j] , 

        for j a power of two, j < self.INT_FIELDS, if the bit
        of valence j in i is equal to 0.

        Each result is reduced modulo p.
        """
        bsh = j << self.LOG_FIELD_BITS
        msk = self.H_MASK[j] 
        t = self.vars.temp() 
        if 2 * bsh  == self.INT_BITS:
            # We may save a bit in this case
            all_mask = self.ALL1_MASK
            s = t.assign((var << bsh) | (var >> bsh))
            self.n_operations += 3
        else:
            s = t.assign(((var << bsh) & msk) | ((var & msk) >> bsh))
            self.n_operations += 5
        s += var.assign((var ^ msk) + t)
        self.n_operations += 2
        self.matrix_code += s
        self.n_code_lines += 2 
        self.reduce_butterfly(var)


    def external_butterfly(self, v1, v2):
        """Auxilary function for method butterfly_op().
   
        Perform butterfly operation on the two C variables v1
        and v2.  For all components v1[i], v2[i] of of v1, v2
        (stored in v1 and v2) we put

            v1[i], v2[i] = v1[i] + v2[i], v1[i] - v2[i] , 

        Each result is reduced modulo p.
        """
        mask = self.H_MASK_EXTERN
        t = self.vars.temp()
        s = t.assign( v1 + (v2 ^ mask) )
        s += v1.assign(v1 + v2)
        self.matrix_code += s
        self.n_code_lines += 2 
        self.n_operations += 3
        self.reduce_butterfly(t, v2)
        self.reduce_butterfly(v1)

    def external_butterfly_all(self, j):
        """Auxilary function for method butterfly_op().
   
        Performs a complete butterfly operation, where each single
        operation comprises two C variables.

        For all components v[i] of v  we put

            v[i], v[i+j] = v[i] + v[j], v[i] - v[j] , 

        for j a power of two, j >= self.INT_FIELDS, if the bit
        oc valence j in i is qualt to 0.
        """
        if self.verbose: 
            print("ext butterfly", j, self.INT_FIELDS)
        assert j >= self.INT_FIELDS
        y = j >> self.LOG_INT_FIELDS
        for i in  range(len(self.vars)):
            if i & y == 0:
                self.external_butterfly(self.vars[i], self.vars[i + y])

    def hadamard_op(self, operations = -1):
         """Multiply vector v by Hadamard matrix H.

         Internal operation:

         For sh = 0, ..., self.LOG_VLEN - 1 we do the following:
         Let y = 1 << sh. Then we put

           v[i], v[i+y]  = (v[i] + v[y]) % p , (v[i] - v[i+y]) % p

         for all i with i & y == 0.

         Remerk on internal operation:

         In cases p = 3, 15, 255, perform method expand_hadamard()
         before doing any butterfly operations. After that 
         expansion, sh has to run through  1,...,self.LOG_VLEN 
         instead. This expansion is reversed after the operation,
         so that the user need not bother about this special case.  
         """
         operations &= (1 << self.LOG_VLEN) - 1
         pre_expanded = self.expanded
         if not self.expanded:
             self.expand_hadamard(self.LOG_VLEN)
         expand_ = (operations & 1) << self.expanded
         operations &=  ~1 & ~(1 << self.expanded)
         operations |=  expand_
         for lg_i in range(self.LOG_VLEN + 1):
             i = 1 << lg_i
             if i & operations == 0:
                 continue
             cmt = "Butterfly: v[i], v[i+%d] = v[i]+v[i+%d], v[i]-v[i+%d]" 
             self.comment(cmt  % (i, i, i))
             if i < self.INT_FIELDS:
                 for var in self.vars:
                     self.internal_butterfly(var, i)     
             else:    
                 self.external_butterfly_all(i)
             self.comment_vector()
         if not pre_expanded:
            self.compress_hadamard()


    def expand_hadamard_intern(self, v1, position):
        """Auxiliary method method for method expand_hadamard()

        It operats on a variables v1 type uint_mmv_t.
        It zeros the components v1[2*i + 1] in vector v1 and it
        copies the old component v1[2*i + 1] to component 
        v1[2*i + 2**position].
        Components v1[2*i + 1] of v1 are set to zero.
        """
        assert 1 <= position < self.LOG_INT_FIELDS
        SH = (self.FIELD_BITS << position) - self.FIELD_BITS
        LO_MASK = self.B0_MASK[1] & self.B0_MASK[1 << position] 
        HI_MASK = LO_MASK << self.FIELD_BITS
        s = v1.assign((v1 & LO_MASK) | ((v1 & HI_MASK) << SH))
        self.matrix_code += s
        self.n_code_lines += 1
        self.n_operations += 4

    def expand_hadamard_intern_all(self, position):
        """Auxiliary method method for method expand_hadamard()

        It executes method expand_hadamard() in case 
        1 <= position < self.LOG_INT_FIELDS.
        """
        for var in self.vars:
            self.expand_hadamard_intern(var, position)


    def expand_hadamard_extern(self, v1, v2):
        """Auxiliary method method for method expand_hadamard()

        It operats on two variables v1, v2 of type uint_mmv_t.
        It zeros the components v1[2*i + 1] in vector v1 and it
        copies the old component v1[2*i + 1] to component v2[2*i].
        Components v1[2*i + 1] of v1 are set to zero.
        """
        mask = self.H_MASK_EXTERN
        s = v2.assign( (v1 >> self.FIELD_BITS) & mask )
        s += v1.assign( v1 & mask )
        self.n_code_lines += 2
        self.n_operations += 3
        self.matrix_code += s


    def expand_hadamard_extern_all(self, position):
        """Auxiliary method method for method expand_hadamard()

        It executes method expand_hadamard() in case 
        position >= self.LOG_INT_FIELDS.
        """
        assert self.LOG_INT_FIELDS <= position <= self.LOG_VLEN
        if position == self.LOG_VLEN:
            self.vlen <<= 1
            self.vars.resize(self.vlen)
        s = ""
        j = 1 << (position - self.LOG_INT_FIELDS)
        mask = self.H_MASK_EXTERN
        for i in range(self.vlen):
            if i & j == 0:
                v1, v2 = self.vars[i], self.vars[i + j]
                self.expand_hadamard_extern(v1, v2)
                

    def expand_hadamard(self, position = -1):
        """This is a very peculiar function for p = 2**2**j-1. 

        Then each component v[i] of v is 2**j bits long. E.g. in case 
        p = 3 we have j = 1 and v[i] is 2 bits long. 
        For a butterfly operation we need one more bit to store the carry
        bit after an addition or subtraction of two components v[i], v[j].
        Thus a butterfly operation cannot be done immediately after 
        loading a vector v into C variables.
                
        This function splits one C variable into two parts, so that 2*2**j
        bits are available for each component v[i], and hence butterfly 
        operations can be done. Use method compress_characteristic_p()
        to undo the effect of this function before storing the vector v to
        to an array.
        
        A slightly nasty side effect is that this function changes the
        order of the components of v as follows:

        Let l be the length of vector v. After calling this method 
        the components of v are ordered a follows:

           0, 2, ..., l-2,   1, 3, ..., l-1 .

        At present this method is used by method hadamard_op()
        only. Method hadamard_op() takes care of this change
        and it restores the vector to its original form.
        So at present the user need not bother about this case.

        By default, entry  2*i+1  is moved to position
        2*i + 2**self.LOG_VLEN.  If the optional argument
        'position' satisfies 1 <= position < self.LOG_VLEN
        then entry  2*i+1  is moved to position
        2*i + 2**position  instead, overwriting the enty at
        thas position. This may be useful in cases where 
        the entries at those positions are not nneded.
        """ 
        if self.NO_CARRY and not self.expanded:
            if not 0 < position < self.LOG_VLEN: 
                 position = self.LOG_VLEN
            if self.verbose:   
               print("expand vector for Hadamard operation")          
            self.comment(
            """Expansion for Hadamard operation:
There is no space for a carry bit between bit fields. So 
we move bit field 2*i + 1  to bit field 2*i + %d.""" % (
                  1 << position)
            )
            if position < self.LOG_INT_FIELDS:             
                self.expand_hadamard_intern_all(position)
            else:
                self.expand_hadamard_extern_all(position)
            self.comment_vector()
            self.expanded = position
        

    def do_compress_hadamard_intern(self, position):
        """Auxiliary method method for method compress_hadamard()

        It executes method compress_hadamard() in case 
        position < self.LOG_INT_FIELDS.
        """
        assert 1 <= position < self.LOG_INT_FIELDS
        TOTAL_BITS = self.FIELD_BITS << position
        SH = TOTAL_BITS - self.FIELD_BITS
        LO_MASK = self.B0_MASK[1] & self.B0_MASK[1 << position] 
        HI_MASK = LO_MASK << self.FIELD_BITS
        for var in self.vars:
            s = var.assign( (var & LO_MASK) | ((var >> SH) & HI_MASK) )
            self.matrix_code += s
            self.n_code_lines += 1
            self.n_operations += 4
 

    def compress_hadamard_extern(self, v1, v2):
        """Auxiliary method method for method compress_hadamard()

        It operates on two variables v1, v2 of type uint_mmv_t. It
        revereses the effect of method expand_hadamard_extern().
        """
        s = v1.assign_xor( v2 << self.FIELD_BITS )
        self.n_code_lines += 1
        self.n_operations += 2
        self.matrix_code += s


    def do_compress_hadamard_extern(self, position):
        """Auxiliary method method for method compress_hadamard()

        It executes method compress_hadamard() in case 
        position >= self.LOG_INT_FIELDS.
        """
        assert self.LOG_INT_FIELDS <= position <= self.LOG_VLEN
        s = ""
        j = 1 << (position - self.LOG_INT_FIELDS)
        mask = self.H_MASK_EXTERN
        for i in range(self.vlen):
            if i & j == 0:
                v1, v2 = self.vars[i], self.vars[i + j]
                self.compress_hadamard_extern(v1, v2)
        if position == self.LOG_VLEN:
            self.vlen >>= 1
            self.vars.resize(self.vlen)
        self.expanded = 0 


    def compress_hadamard(self):
        """Reverse the effect of method expand_hadamard().""" 
        if self.expanded: 
            if self.verbose:   
                print("compress vector after Hadamard operation")         
            self.comment("Reverse expansion for Hadamard operation")
            if self.expanded < self.LOG_INT_FIELDS:             
                self.do_compress_hadamard_intern(self.expanded)
            else:
                self.do_compress_hadamard_extern(self.expanded)
            self.expanded = 0
            self.comment_vector()
 

    def mul_var_pwr2(self, var, l):
        """Multiply all components of variable 'var' by scalar 2**l.

        Here l may be any integer.
        Modulo p = 2**k-1, this is simply a bit rotation.
        """
        ls = l % self.P_BITS
        if ls != 0:
            hs = self.P_BITS - ls
            m_l = self.smask(range(hs))
            m_h = self.smask(range(hs, self.P_BITS))
            s = var.assign(((var & m_l) << ls) | ((var & m_h) >> hs))
            self.matrix_code += s
            self.n_code_lines += 1
            self.n_operations += 5

      

    def mul_pwr2(self, l):
        """Multiply vector v by the scalar 2**l.

        Here l may be any integer.
        Modulo p = 2**k-1, this is simply a bit rotation.
        """
        if l % self.P_BITS == 0:
            s = "Multiplication by 2**%d is trivial mod %d"
            self.comment(s % (l, self.P))
            return
        s = "Multiply vector by scalar 2**%d mod %d"
        self.comment(s % (l, self.P))
        for w in self.vars:
            self.mul_var_pwr2(w, l)
            


    def load_vector_direct(self, array_name):
        """Load vector v form an (external) array A of integers.

        'array_name' must a variable of type uint_mmv_t 
        referring to array A.
        """
        self.comment("Loading vector v from array %s" % array_name)
        s = ""
        for i,v in enumerate(self.vars):
            s += "%s = %s[%d];\n" % (v, array_name, i) 
            self.n_code_lines += 1
        self.matrix_code += s
        self.comment_vector()  

    def store_vector_direct(self, array_name):
        """Store a vector v to an (external) array A of integers.

        'array_name' must a variable of type uint_mmv_t 
        referring to array A.
        """ 
        self.comment("Storing vector v to array %s" % array_name)
        s = ""
        for i,v in enumerate(self.vars):
            s += "%s[%d] = %s;\n" % (array_name, i, v) 
            self.n_code_lines += 1
        self.matrix_code += s


    def complement_variable(self, var):
        """Generate code: 'var = ~var'; for variable 'var'"""
        self.code("%s = ~(%s);\n" % (var, var))
        self.n_code_lines += 1
        self.n_operations += 1

    def code(self, string):
        self.matrix_code += string;



    def generate(self):
        if self.verbose:
            print("Generate C code using class HadamardMatrixCode")
        """Return C code generated so far as a string"""
        s = """
// This is an automatically generated matrix operation, do not change!
{
"""     
        s += self.vars.declare() + "\n"
        s += self.matrix_code 
        s +="""}
// End of automatically generated matrix operation.
 
"""   
        self.reset_vars()
        return s

