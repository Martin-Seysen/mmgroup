

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import re
import numpy as np
import warnings
from functools import partial

from mmgroup.structures.abstract_mm_group import AbstractMMGroupWord
from mmgroup.structures.abstract_mm_group import AbstractMMGroup
from mmgroup.structures.parse_atoms import AtomDict, ihex     
from mmgroup.structures.construct_mm import iter_mm       

from mmgroup.structures.ploop import Cocode, PLoop
from mmgroup.structures.autpl import StdAutPlGroup, AutPL ,autpl_from_obj
from mmgroup.structures.construct_mm import iter_mm       
from mmgroup.structures.construct_mm import load_group_name     

from mmgroup.mat24 import ploop_theta
from mmgroup.generators import gen_leech2_type
from mmgroup.clifford12 import xsp2co1_elem_to_qs_i, xsp2co1_elem_to_qs 
from mmgroup.clifford12 import xsp2co1_qs_to_elem_i 
from mmgroup.clifford12 import xsp2co1_chain_short_3 
from mmgroup.clifford12 import xsp2co1_neg_elem
from mmgroup.clifford12 import error_string, chk_qstate12
from mmgroup.clifford12 import xsp2co1_unit_elem
from mmgroup.clifford12 import xsp2co1_mul_elem, xsp2co1_inv_elem
from mmgroup.clifford12 import xsp2co1_copy_elem
from mmgroup.clifford12 import xsp2co1_reduce_elem
from mmgroup.clifford12 import xsp2co1_elem_to_leech_op
from mmgroup.clifford12 import xsp2co1_set_elem_word 
from mmgroup.clifford12 import xsp2co1_mul_elem_word 
from mmgroup.clifford12 import xsp2co1_xspecial_vector
from mmgroup.clifford12 import xsp2co1_xspecial_conjugate
from mmgroup.clifford12 import xsp2co1_elem_xspecial
from mmgroup.clifford12 import xsp2co1_elem_to_word
from mmgroup.clifford12 import xsp2co1_involution_invariants
from mmgroup.clifford12 import xsp2co1_involution_orthogonal
from mmgroup.clifford12 import xsp2co1_conjugate_elem
from mmgroup.clifford12 import xsp2co1_elem_conjugate_2B_involution

from mmgroup.structures.qs_matrix import QStateMatrix


FORMAT_REDUCED = True



###########################################################################
# Word class for the group G_{x0}
###########################################################################


MM = None

def import_MM():
    global MM
    if MM is not None:
        return MM
    from mmgroup.mm_group import MM
    return MM


class Xsp2_Co1(AbstractMMGroupWord):
    """Model an element the subgroup :math:`G_{x0}` of the Monster

    See class ``Xsp2_Co1`` for the definition of that group.   

    The constructor of this class returns the neutral element of
    the group :math:`G_{x0}`. 

    :var group:
        This attribute contains the group to which the element belongs.
        That group is an instance of class Xsp2_Co1.

    .. warning::
       The constructor of this class is not for public use! You
       may call an instance ``M`` of class Xsp2_Co1 for
       constructing elements of the instance ``M`` of the monster 
       group.
  
    """
    MIN_LEN = 16
    __slots__ =  "_data"
    def __init__(self,  tag = None, atom = None, *args, **kwds):
        atoms = iter_mm(self.group, tag, atom, in_G_x0 = True)
        a_atoms = np.fromiter(atoms, dtype = np.uint32) 
        self._data = np.zeros(26, dtype = np.uint64)
        xsp2co1_set_elem_word(self._data, a_atoms, len(a_atoms))
         
    @property
    def data(self):
        return list(map(int, self._data))

    @property
    def short3(self):
        return int(self._data[0])

    @property
    def qs_t(self):
        return QStateMatrix(xsp2co1_elem_to_qs_i(self._data))
        
    @property
    def qs(self):
        return QStateMatrix(xsp2co1_elem_to_qs(self._data))
        
    @property
    def leech_op(self):
        a = np.zeros(576, dtype = np.int8)
        xsp2co1_elem_to_leech_op(self._data, a) 
        return a.reshape((24,24))        

        
    def order(self, max_order = 119):
        """Return the order of the element of the monster group


        If the argument ``max_order`` is present then the order of the 
        element is checked up to (and including) ``max_order`` only.  
        Then the function returns ``0`` if the order is greater than 
        ``max_order``. By default, the function returns the exact 
        order of the element.
        """
        o = self.qs.order(max_order)
        if o & 1 == 0:
            o = o >> 1
        unit, pw = self.group(), self**o
        if pw == unit:
            return o
        for i in range(2):
            o, pw = 2*o, pw * pw
            if pw == unit:
                return o
        err = "Order of QStateMatrix object not found" 
        raise ValueError(err)


    def as_xsp(self):
        return chk_qstate12(xsp2co1_xspecial_vector(self._data))

    as_Q_x0_atom = as_xsp

    def type_Q_x0(self):
        """Return type of element if it is in the subgroup :math:`Q_{x0}`

        If the element is in the subgroup :math:`Q_{x0}` of the monster
        then the function returns the type of the vector in the Leech 
        lattice modulo 2 corresponding to this element. That type is
        0, 2, 3, or 4.

        The function raises ValueError if the element is not
        in the subgroup :math:`Q_{x0}`. 
        """
        v = self.as_Q_x0_atom()
        return gen_leech2_type(v) >> 4

    def xsp_conjugate(self, v, sign = True):
        v = np.array(v, dtype = np.uint64, copy=True)
        shape = v.shape
        assert len(shape) <= 1
        v = np.ravel(v)
        chk_qstate12(
            xsp2co1_xspecial_conjugate(self._data, len(v), v, sign))
        if len(shape):
            return list(map(int,v))
        else:
            return int(v[0])

    def mul_data(self, data):
        a_atoms = np.array(data, dtype = np.uint32)
        xsp2co1_mul_elem_word(self._data, a_atoms, len(a_atoms))
        return self


    @property
    def mmdata(self):
        a = np.zeros(10, dtype = np.uint32)
        length = chk_qstate12(xsp2co1_elem_to_word(self._data, a))
        return a[:length]

    def _involution_invariants(self):
        """Wrapper for C function  xsp2co1_involution_invariants"""
        invar = np.zeros(12, dtype = np.uint64)
        xsp2co1_involution_invariants(self._data, invar)
        v1 = xsp2co1_involution_orthogonal(invar, 1)
        v0 = xsp2co1_involution_orthogonal(invar, 0)
        return invar, v1, v0


    def conjugate_2B_involution(self, mmgroup = None):
        """Find an element conjugating an involution into the centre

        If the element :math:`g` given by ``self`` is a 2B involution 
        in  the monster then the method computes an element :math:`h` 
        of the monster   with  :math:`h^{-1} g h = z`, where  
        :math:`z` is the central involution in the 
        subgroup :math:`G_{x0}` of the  monster.

        The function returns  :math:`h` as an element of the instance
        ``MM`` of class ``MMGroup``. It raises ``ValueError``
        if :math:`g` is not a 2B involution in the monster. 

        This is a wrapper for the C 
        function ``xsp2co1_elem_conjugate_2B_involution``.

        If parameter ``mmgroup`` is an instance of class ``MMGroup``
        then the function returns :math:`h` as an element of the 
        group given by ``mmgroup``.
        """
        if mmgroup is None:
             mmgroup = MM if MM is not None else import_MM()            
        a = np.zeros(14, dtype = np.uint32)
        len_a = xsp2co1_elem_conjugate_2B_involution(self._data, a)
        chk_qstate12(len_a)
        return mmgroup('a', a[:len_a])



###########################################################################
# The class representing the group G_x0
###########################################################################



def cocode_to_xsp2co1(g, c):
    res =  g.word_type(group = g)
    chk_qstate12(xsp2co1_elem_xspecial(res._data, c.value & 0xfff))
    return res

def ploop_to_xsp2co1(g, pl):
    res =  g.word_type(group = g)
    value = (c.value & 0x1fff)
    value = (value << 12) ^ ploop_theta(value)
    chk_qstate12(xsp2co1_elem_xspecial(res._data, value))
    return res


def autpl_to_xsp2co1(g, aut):
    res =  g.word_type(group = g)
    a = np.zeros(2, dtype = uint32)
    a[0] = 0x10000000 + (aut._cocode & 0xfff)   # tag 'd'
    a[1] = 0x20000000 + aut._perm_num           # tag 'p'
    chk_qstate12(xsp2co1_set_elem_word(res._data, a, 2))
    return res

def mmgroup_to_xsp2co1(g, mm):
    res =  g.word_type(group = g)
    chk_qstate12(xsp2co1_set_elem_word(res._data, mm._data, mm.length))
    return res


def xsp2co1_to_mm(mmgroup, xsp):
    g = mmgroup()
    g._extend(10)
    g.length = chk_qstate12(xsp2co1_elem_to_word(xsp._data, g._data))
    g.reduced = 0 
    return g



class Xsp2_Co1_Group(AbstractMMGroup):
    r"""Model the subgroup :math:`G_{x0}` of the Monster
    
    The group :math:`G_{x0}` is the subgroup of the Monster group
    genertated by the elements with tags ``p, d, x, y, l`` in the  
    class |MMGroup| representing the Monster group.  :math:`G_{x0}`
    has structure :math:`2**(1+24).\mbox{Co}_1` in ATLAS
    notation, see :cite:`Con85`, :cite:`Atlas`.
       

    
    :param \*data:

      A variable number of arguments; each argument describes an
      element of the subgroup ``2**(1+24).Co_1`` of the monster.
      These elements are multiplied.  
        

    Depending on its type each parameter in **\*data** is  
    interpreted as follows:

    .. table:: Legal types for constructor of class ``Xsp2_Co1_Group``
      :widths: 25 75

      ===================== ==================================================
      type                  Evaluates to
      ===================== ==================================================
      tuple (``tag, data``) Create an element of the monster group as 
                            described in class |MMGroup|. Here legal tags
                            are: ``'p', 'd', 'x', 'y', 'z', 'l'``.  The
                            resulting element must lie in the subgroup
                            :math:`G_{x0}` of the monster.
                            
      class |MMGroup|       Create an element of the Monster group. That                
                            element must lie in the subgroup
                            :math:`G_{x0}` of the monster.

      class |AutPL|         Create an element of the subgroup of
                            :math:`G_{x0}` given by class |AutPL|.
                           
      class ``Xsp2_Co1``    Create a copy of an element of this class.
                           
      pair (``qstate, x``)  Deprecated and not implemented!!!!!
                            This kind of construction is for testing and
                            not recommended for public use. Here        
                            ``qstate`` must be an instance of class 
                            ``mmgroup.structures.qs_matrix.QStateMatrix``
                            and ``x`` must be an integer representing a
                            short Leech lattice vector modulo 3.
      ===================== ==================================================

        
    :raise:
        * TypeError if ``type(data)`` is not as expected.
        * ValueError if ``data`` cannot be converted to an
          instance of class  `` Xsp2_Co1``.
    
    
    """
    __instance = None
    __slots__ = "data"
    STD_V3  = 0x8000004
    word_type = Xsp2_Co1
    is_mmgroup = True


    def __new__(cls):
        if Xsp2_Co1_Group.__instance is None:
             Xsp2_Co1_Group.__instance = AbstractMMGroup.__new__(cls)
        return Xsp2_Co1_Group.__instance


    def __init__(self):
        """ TODO: Yet to be documented     


        """
        super(Xsp2_Co1_Group, self).__init__()
        self.atom_parser = AtomDict(self.atom)


    def atom(self, tag = None, i = "r"):
        return self.word_type(tag, i)

    def _imul(self, g1, g2):
        chk_qstate12(xsp2co1_mul_elem(g1._data, g2._data, g1._data))
        return g1

    def _invert(self, g1):
        w = self.word_type()
        chk_qstate12(xsp2co1_inv_elem(g1._data, w._data))
        return w

    def copy_word(self, g1):
        w = self.word_type()
        xsp2co1_copy_elem(g1._data, w._data)
        return w

    def reduce(self, g1, copy = False):
        if copy:
            g1 = self.copy_word(g1)
        chk_qstate12(xsp2co1_reduce_elem(g1._data))
        return self
       
    def _equal_words(self, g1, g2):
        chk_qstate12(xsp2co1_reduce_elem(g1._data))
        chk_qstate12(xsp2co1_reduce_elem(g2._data))
        return (g1._data == g2._data).all()

    def _embed_number(self, n):
        w = self.word_type()
        if (n == -1):
            xsp2co1_neg_elem(w._data)
            n = 1
        if n == 1:
            return w
        raise TypeError("Cannot convert a number to a group element")
        
    def from_qs(self, qs, x):  
        w = self.word_type()
        w0 =  xsp2co1_qs_to_elem_i (qs, x)
        for i in range(26):
             w._data[i] =  w0[i]
        return w             

    def str_word(self, v1):
        return "Xsp2_Co1 " + str_xsp2_co1(v1._data)
 
    def from_xsp(self, x):
        w = self.word_type()
        chk_qstate12(xsp2co1_elem_xspecial(w._data, x))
        return w

    def from_data(self, data):
        """Create a group element from an array of generators

        Internally, an element of group is represented
        as an array of unsigned 32-bit integers, where each entry
        of the array describes a generator. See section
        :ref:`header-mmgroup-generators-label` for details.
 
        This function creates an element of the group from
        such an array of integers.

        :param data: An array-like object representing a 
                     word of generators of the monster group

        :return: An element of this instance of the group 
        :rtype:  an instance of class 
                 mmgroup.structures.xsp2_co1.Xsp2_Co1_Word

        """
        return self.word_type('a', data)




_dict_pm3 = {0: '0', 1:'+', 0x1000000:'-', 0x1000001:'0'}
def str_leech3(x):
    x = int(x)
    lst = [_dict_pm3[(x >> i) & 0x1000001] for i in range(24)]
    return "(" + "".join(lst) + ")"

def str_xsp2_co1(data, factor = 1, t = False):
    qs0 = xsp2co1_elem_to_qs_i(data) if t else xsp2co1_elem_to_qs(data)
    qs = QStateMatrix(qs0) / factor
    return str_leech3(data[0]) + " (x) " + str(qs)


try:
    from mmgroup.clifford12 import xsp2co1_error_pool
    def get_error_pool(length):
        assert length > 0
        a = np.zeros(length, dtype = np.uint64)
        length = xsp2co1_error_pool(a, length)
        return list(map(int, a[:length])) 
except:
    def get_error_pool(length):        
        return []    





StdXsp2_Co1_Group = Xsp2_Co1_Group()
Xsp2_Co1.group = StdXsp2_Co1_Group
load_group_name(StdXsp2_Co1_Group, "MX")



   