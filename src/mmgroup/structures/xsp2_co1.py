

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import re
import numpy as np
import warnings

from mmgroup.structures.auto_group import AbstractGroupWord
from mmgroup.structures.auto_group import AbstractGroup
from mmgroup.structures.parse_atoms import AtomDict, ihex     

from mmgroup.structures.ploop import Cocode, PLoop
from mmgroup.structures.autpl import StdAutPlGroup, AutPL ,autpl_from_obj



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

from mmgroup.structures.qs_matrix import QStateMatrix

from mmgroup.mm_group import gen_atom, MM

FORMAT_REDUCED = True



###########################################################################
# Word class for the group G_{x0}
###########################################################################



class Xsp2_Co1_Word(AbstractGroupWord):
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
    def __init__(self, atoms = [], **kwds):
        self.group = kwds['group']
        self._data = np.zeros(26, dtype = np.uint64)
        a_atoms = np.array(atoms, dtype = np.uint32)
        xsp2co1_set_elem_word(self._data, a_atoms, len(a_atoms))
        """
        if len(atoms):
            chk_qstate12(xsp2co1_mul_elem_atom(self._data, atoms[0], 1))
            for v in atoms[1:]:
                chk_qstate12(xsp2co1_mul_elem_atom(self._data, v, 0))
        else:
            xsp2co1_unit_elem(self._data)
        """

         
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

    def xsp_conjugate(self, v):
        v = np.array(v, dtype = np.uint64, copy=True)
        shape = v.shape
        assert len(shape) <= 1
        v = np.ravel(v)
        chk_qstate12(xsp2co1_xspecial_conjugate(self._data, len(v), v))
        if len(shape):
            return list(map(int,v))
        else:
            return int(v[0])

    def mul_data(self, data):
        a_atoms = np.array(data, dtype = np.uint32)
        xsp2co1_mul_elem_word(self._data, a_atoms, len(a_atoms))
        return self



###########################################################################
# The class representing the group G_x0
###########################################################################



def cocode_to_xsp2co1(g, c):
    res =  g.word_type(group = g)
    chk_qstate12(xsp2co1_op_delta_pi(res._data, c.cocode, 0))
    return res

def autpl_to_xsp2co1(g, aut):
    res =  g.word_type(group = g)
    chk_qstate12(xsp2co1_op_delta_pi(res._data, aut.cocode, aut.perm_num))
    return res

def mmgroup_to_xsp2co1(g, mm):
    raise NotImplementedError


class Xsp2_Co1_Group(AbstractGroup):
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

    .. table:: Legal types for constructor of class ``AutPL``
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
    __slots__ = "data"
    STD_V3  = 0x8000004
    word_type = Xsp2_Co1_Word
    tags, formats = " dpxyl", [None, ihex, str, ihex, ihex, str]
    atom_parser = {}               # see method parse()
    conversions = {
        Cocode: cocode_to_xsp2co1,
        AutPL: autpl_to_xsp2co1,
    }
    FRAME = re.compile(r"^M?\<(.+)\>$") # see method parse()
    STR_FORMAT = r"M<%s>"

    def __init__(self):
        """ TODO: Yet to be documented     


        """
        super(Xsp2_Co1_Group, self).__init__()
        self.atom_parser = AtomDict(self.atom)
        self.set_preimage(StdAutPlGroup,  tuple)


    def atom(self, tag = None, i = "r"):
        return self.word_type(gen_atom(tag, i), group = self)

    def _imul(self, g1, g2):
        chk_qstate12(xsp2co1_mul_elem(g1._data, g2._data, g1._data))
        return g1

    def _invert(self, g1):
        w = self.word_type([], group=self)
        chk_qstate12(xsp2co1_inv_elem(g1._data, w._data))
        return w

    def copy_word(self, g1):
        w = self.word_type([], group=self)
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
        w = self.word_type([], group=self)
        if (n == -1):
            xsp2co1_neg_elem(w._data)
            n = 1
        if n == 1:
            return w
        raise TypeError("Cannot convert a number to a group element")
        
    def from_qs(self, qs, x):  
        w = self.word_type([], group=self)
        w0 =  xsp2co1_qs_to_elem_i (qs, x)
        for i in range(26):
             w._data[i] =  w0[i]
        return w             

    def str_word(self, v1):
        return "Xsp2_Co1 " + str_xsp2_co1(v1._data)
 
    def from_xsp(self, x):
        w = self.word_type([], group=self)
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
        return self.word_type(data, group = self)



Xsp2_Co1 = Xsp2_Co1_Group()
Xsp2_Co1.set_preimage(MM, tuple)
MM.set_preimage(Xsp2_Co1, tuple)



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



   