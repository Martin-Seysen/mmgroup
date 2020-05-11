"""Access to the 196884-dimensional representation of the monster"""

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import sys
import os
import numpy as np
from numbers import Integral
from random import randint
import warnings
import time
from importlib import import_module






from mmgroup.structures.abstract_mm_rep_space import AbstractMmRepVector
from mmgroup.structures.abstract_mm_rep_space import AbstractMmRepSpace
from mmgroup.mm_group  import MMGroup, MMGroupWord


from mmgroup.mm import mm_vector, mm_aux_random_mmv
from mmgroup.mm import mm_aux_zero_mmv, mm_aux_reduce_mmv
from mmgroup.mm import mm_aux_mmv_to_sparse, mm_aux_mmv_add_sparse
from mmgroup.mm import mm_aux_mmv_extract_sparse
from mmgroup.mm import mm_aux_mmv_set_sparse, mm_aux_mmv_add_sparse
from mmgroup.mm import mm_aux_bytes_to_mmv, mm_aux_mmv_to_bytes
from mmgroup.mm import mm_aux_get_mmv, mm_aux_put_mmv
from mmgroup.mm import mm_rng_make_seed
from mmgroup.mm import INT_BITS, PROTECT_OVERFLOW
from mmgroup.mm import mm_aux_check_mmv
from mmgroup.mm import mm_aux_get_mmv1
from mmgroup.mm import mm_sparse_purge


uint_mmv = np.uint32 if INT_BITS == 32 else np.uint64
standard_seed = mm_rng_make_seed()
standard_mm_group = MMGroup()



######################################################################
# Importing a C wrapper for a specific characteristic 'p'
######################################################################

mm_op = {}
all_characteristics_found = None

def mm_wrapper(p):
    """Return the module dealing with 'characteristic' p"""
    try:
        return mm_op[p]
    except KeyError:
        mm_op[p] = import_module('mmgroup.mm%d' % p)
        return mm_op[p] 

def characteristics():
    """Return list of all 'characteristics' p supported"""
    global all_characteristics_found
    if all_characteristics_found is None:
        for k in range(2,8):
            p = (1 << k) - 1
            if not p in mm_op:
                try:
                    mm_op[p] = import_module('mmgroup.mm%d' % p)
                except ModuleNotFoundError:
                    pass
        all_characteristics_found = sorted(mm_op.keys())
    return all_characteristics_found

######################################################################
# Modelling a vector of the 196884-dimensional rep of the monster
######################################################################

class MMSpaceVector(AbstractMmRepVector):
    def __init__(self, space):
        self.space = space
        self.data = mm_vector(self.space.p)

    def check(self, comment = ""):
        """Check if the vector is correct

        Raise ValueError if the vector is errorneous.
        """
        self.space.check(self, comment = "")

       
    def mul_exp(self, g, e, break_g = False):
        """Multiply the vector with g**e inplace

        The vector is changed and the changed vector is returned.
        
        Afterwards, the vector contains an attribute 'last_timing'
        containing the run time of this operaiton in seconds,
        """
        return self.space.vector_mul_exp(self, g, e, break_g)
 
        

######################################################################
# class MMSpace
######################################################################




class MMSpace(AbstractMmRepSpace):
    """Models the sparse representation 198884x of the monster group. 

    YET TO BE DOCUMENTED !!!

    """
    vector_type = MMSpaceVector

    check_errors = {
        -1: "Bad input value p",
        -2: "A one bit outside a field has been found",
        -3: "A subfield has an illegal nonzero entry at index >= 24",
        -4: "Illegal nonzero diagonal entry", 
        -5: "Symmetric part of vector is not symmetric",
    }

    extend_modulus = {
        3:3, 7:7, 15:15, 31:31, 63:63, 127:127, 255:255,
        5:15, 9:63, 17:255,
    }

    tag_indices = { # tag: (start, nrows, ncolumns, stride)
       "A": (     0,   24, 24, 32),
       "B": (   768,   24, 24, 32),
       "C": (  1536,   24, 24, 32),
       "T": (  2304,  759, 64, 32),
       "X": ( 50880, 2048, 24, 32),
       "Z": (116416, 2048, 24, 32),
       "Y": (181952, 2048, 24, 32),
    }

    def __init__(self, p, group = None):
        """Create a 196884-dimensional representation of the monster

        All calculations are done modulo the odd number p
        """
        try: 
            self.p = self.extend_modulus[p]
            self.p_original = p
        except KeyError:
            raise KeyError("Modulus %s not supported for monster group" % p)
        if group is None:
            group = standard_mm_group 
        assert isinstance(group, MMGroup) 
        super(MMSpace, self).__init__(self.p, group)
        mm = mm_wrapper(p)  # fetch mm_op[p]
        #self.mm = mm_wrapper(self.p)
                            # if we do this, we cannot pickle vectors
        self.MMV_INTS = mm_op[self.p].MMV_INTS 
        self.op_vector_add = mm_op[self.p].op_vector_add
        self.op_scalar_mul = mm_op[self.p].op_scalar_mul
        self.op_word = mm_op[self.p].op_word
        self.op_compare = mm_op[self.p].op_compare
        del mm

    @property
    def mm(self):
        """Return module object mmgroup.mm<p> for characteristic p"""
        return mm_op[self.p]


    #######################################################################
    # Creating vectors 
    #######################################################################

    def zero(self):
        """Create zero vector with minimal overhead"""
        return MMSpaceVector(self)

    def copy_vector(self, v1):
        assert v1.space == self
        v = MMSpaceVector(self)
        np.copyto(v.data, v1.data)
        return v

    def rand_uniform(self, seed = None):
        v = MMSpaceVector(self)
        if seed is None:
            seed = standard_seed
        mm_aux_random_mmv(self.p, seed, v.data) 
        return v

    #######################################################################
    # Obtaining and setting components via sparse vectors
    #######################################################################

    def getitems_sparse(self, v1, a_indices):
        """Get items from vector v1

        Here we assert that v1 is a vector of this vector space and
        that 'a_indices' is a one-dimensional numpy array of type
        np.uint32, containing the coordinates to be read from v1.
 
        The function must add the corresponding coordinate to each
        entry of the array 'sparse_items'. All coordinates must be
        nonnegative and < 256. 

        A zero entry in the array 'a_indices' is ignored.
        """
        assert v1.space == self
        if len(a_indices):
            mm_aux_mmv_extract_sparse(self.p, v1.data, a_indices,
                len(a_indices))
        return a_indices 

    def additems_sparse(self, v, a_indices):
        """Add a vector in sparse representation to vector v.

        This method takes a numpy array 'a_indices' of integers of dtype 
        numpy.uint32 containing the description of a vector v2 in sparse 
        representation. It computes 

             v  =  v + v2 .

        Here vector v is a standard vector in this space.
        """
        if len(a_indices):
            mm_aux_mmv_add_sparse(self.p, a_indices, len(a_indices),
                v.data)
        return v

    def setitems_sparse(self, v, a_indices):
        """Setro selected components of a vector 

        Arguments 'v' and 'a_indices' are as in method getitems_sparse().
        Here the coordinates of vector 'v' described by 'a_indices' are 
        set to the values given in 'a_indices'. 
        The array 'a_indices' is not changed.
        """
        assert v.space == self
        if len(a_indices):
            mm_aux_mmv_set_sparse(self.p, v.data, a_indices, 
                len(a_indices))
        return v


    #######################################################################
    # Conversion from and to to sparse representation 
    #######################################################################

    def as_sparse(self, v1):
        sp = np.zeros(196884, dtype = np.uint32)
        length = mm_aux_mmv_to_sparse(self.p, v1.data, sp)
        return sp[:length]


    #######################################################################
    # Vector operations 
    #######################################################################


    def iadd(self, v1, v2):
        self.op_vector_add(v1.data, v2.data)
        return v1
 
    def imul_scalar(self, v1, a):
        self.op_scalar_mul(a % self.p, v1.data)
        return v1
           
    #######################################################################
    # Group operation 
    #######################################################################

    def imul_group_word(self, v1, g):
        """Return product v1 * g of vector v1 and group word g.

        v1 may be destroyed.

        This method is called for elements v1 of the space
        'self' and for elements g of the group 'self.group' only.
        """
        work =  mm_vector(self.p)
        assert v1.space == self
        g = self.group(g)
        assert isinstance(g, MMGroupWord) 
        self.op_word(v1.data, g._data, g.length, 1, work)
        return v1  


    def vector_mul_exp(self, v1, g, e, break_g = False):
        """Compute product v1 * g**e of vector v1 and group word g.

        Here v1 is a vector in this space, e is an integer, g is a 
        group element, and  v1 is replaced by v1 * g**e.

        This method should be  called for elements v1 of the space
        'self' and for elements g of the group 'self.group' only.

        If break_g is True, each factor g is multiplied with v1
        separately. Otherwise, the expression  g**e  may be
        optimized. This option is mainly for benchmarking.

        After applying this function to vecter v1, the vector
        v1 has an attribute v1.last_timing containing the run
        time of the C part of this operation in seconds.
        """
        work = mm_vector(self.p)
        assert v1.space == self
        assert -1 << 31 < e < 1 << 31
        g = self.group(g)
        assert isinstance(g, MMGroupWord) 
        length = g.length
        if break_g:
             g._extend(length + 1)
             g[length] = 0x70000000
             length += 1
        t_start = time.process_time()
        self.op_word(v1.data, g._data, length, e, work)
        v1.last_timing = time.process_time() - t_start
        return v1  


    #######################################################################
    # Checking equality
    #######################################################################

    def equal_vectors(self, v1, v2):
        """Return True iff vectors v1 and v2 are equal 

        This method is called for elements v1 and v2 of the space
        'self' only.
        """
        if self.p == self.p_original:
            return not self.op_compare(v1.data, v2.data) 
        return self.as_bytes(v1) == self.as_bytes(v2)

    #######################################################################
    # Conversion from and to byte format
    #######################################################################

    def as_bytes(self, v1):
        """Return vector 'self' as a byte array

        The result is a numpy array with dtype = uint8 and
        shape = (196884,).
        """
        b = np.zeros(196884, dtype = np.uint8)
        mm_aux_mmv_to_bytes(self.p, v1.data, b)
        if self.p != self.p_original:
            b %= self.p_original
        return b

    def from_bytes(self, b):
        b = np.array(b, dtype = np.int32)
        if len(b.shape) != 1:
            raise TypeError("Bad shape of byte data vector")
        if len(b) != 196884:
            raise TypeError("Bad length of byte data vector")
        b = np.array(b % self.p, dtype = np.uint8)
        v =self.zero()
        mm_aux_bytes_to_mmv(self.p, b, v.data)
        return v

 
        
    #######################################################################
    #  Checking ad reducing a vector
    #######################################################################

    def check(self, v1, comment = ""):
        """Check the vector 'self'.

        Raise ValueError if an error is found in vector 'self'.
        """
        if comment != "": 
            comment = ", (%s)" % str(comment)
        if len(v1.data) != self.MMV_INTS + 1:
            err = "MM vector has wrong length" + comment
            raise MemoryError(err)   
        if v1.data[-1] != PROTECT_OVERFLOW:
            err = "Buffer overflow in MM vector detected" + comment
            raise MemoryError(err)
        result = mm_aux_check_mmv(self.p, v1.data)
        #print("check result is", result)
        if not result:
            return
        try:
            err = self.check_errors[result] + comment
        except:
            err = "Unknown error %d in MM vector" % result + comment
        print("\n%s!\n" % err)
        raise ValueError("Error in MM vector")
 
    def reduce(self, v1):
        """Convert vector v1 to a unique reduced form"""
        if self.p == self.p_original:
            mm_aux_reduce_mmv(self.p, v1.data)
        else:
            b = np.zeros(196884, dtype = np.uint8)
            mm_aux_mmv_to_bytes(self.p, v1.data, b)
            b %= self.p_original
            mm_aux_bytes_to_mmv(self.p, b, v1.data)
        return v1

     

 