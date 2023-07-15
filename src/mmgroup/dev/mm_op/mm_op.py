r"""We describe the class ``MM_Op``.

Class ``MM_Op`` in module ``mmgroup.dev.mm_op.mm_op`` is a 
table-providing class for the C functions used by python extensions
``mmgroup.mm<p>``. Each of these extensions implements
the operation of monster group modulo a fixed number ``p``. 

Class ``MM_Op`` provides the same functionality as class ``MM_Const``
in module ``mmgroup.dev.mm_basics.mm_basics``. Since modulus ``p``
is fixed in an instance of class ``MM_Op``, all constants provided
by the base class ``MM_Basics`` of ``MM_Op`` and ``MM_Const`` (see 
table *Constants used for generating representations of the monster*)
have fixed values. So they are available as attributes of an
instance of class ``MM_Op`` and they may also be used by the code 
generator directly via string formatting.
 
The string formatting functions ``%{shl:expression,i}`` and 
``%{smask:value, fields, width}`` work as in class ``MM_Op``. In 
addition, parameter ``fields`` defaults to ``-1`` (i.e. the mask 
is set in all bit fields) and ``width`` defaults to the constant 
``FIELD_BITS``. i.e. the number of bits used to store an integer 
modulo ``p``. E.g. for ``p = 7`` we have ``FIELD_BITS = 4``, and 
``%{smask:1}`` evaluates to the ``64``-bit integer constant 
``0x1111111111111111``.

Class ``MM_Op`` also provides the following directives:

``MMV_ROTL src, count, dest``

  Rotate the bit fields of a variable of type ``uint_mmv_t``

  Here ``src`` is an integer of type ``uint_mmv_t`` which is 
  interpreted as an array of bit fields, where each bit field
  stores a number modulo ``p``. Then each bit field is rotated 
  left by ``count`` bits. This means that the numbers in all
  bit fields are multiplied by ``2**count``. ``count`` must be 
  an integer. It is reduced modulo the bit length of ``p``.
  So ``count`` may also be negative. The result of the 
  rotation is written to the variable ``dest`` which should
  also be of type ``uint_mmv_t``. ``dest`` defaults to ``src``.

``MMV_UINT_SPREAD src, dest, value``
  
  Spread bits of an integer ``src`` to the bit fields of ``dest``.
  Here ``src`` may be any integer variable and ``dest`` should
  be a variable of type ``uint_mmv_t``. If bit ``i`` of the
  integer ``src`` is set then the ``i``-th bit field of variable
  ``dest`` is set to ``value``; otherwise it is set to zero.
  ``value`` must be an integer with ``0 <= value <= p``; 
  default is ``p``.

"""

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import sys
import os
import warnings
from numbers import Integral
import re

from mmgroup.generate_c import UserDirective, UserFormat

from mmgroup.dev.mm_basics.mm_basics import INT_BITS, c_hex
from mmgroup.dev.mm_basics.mm_basics import smask, smask_default
from mmgroup.dev.mm_basics.mm_basics import MM_Basics


########################################################################
# Class MM_Op
########################################################################

class MM_Op(MM_Basics):
    r"""Supports basic operations on integers of type ``uint_mmv_t``.

    Here an integer of type ``uint_mmv_t`` is interpreted as an array 
    of bit fields, where each bit field stores an integer modulo ``p``
    for a fixed number ``p``. This class is similar to class 
    ``MM_Const`` in module ``mmgroup.dev.mm_basics.mm_basics``,
    where the modulus ``p`` may be variable.

    Usage of this class is documented in the module documentation.

    :param p: 

       This is the modulus ``3 <= p < 256``. ``p + 1`` must be a 
       power of two.

    :type p: int

    """
    m_psum = re.compile(r"\s*(\w+)\s*([+-].*)?$")

    def __init__(self, **kwds):
        """Create instance for modulus p with p + 1 a power of two"""
        self.p = p = int(kwds.get('p', 3))
        self.directives =  {
            "MMV_ROTL" : UserDirective(self.gen_mmv_rotl, "sis"),
            "MMV_UINT_SPREAD": UserDirective(self.gen_mmv_uint_spread, "ss"),
        }
        self.tables = self.sizes(p).copy()
        self.tables["GENERATE_CODE"] = True
 
    def __getattr__(self, attrib):
        return self.tables[attrib]

    def smask(self, value, fields = -1, width = None):
        """Return smask(value, fields, width)

        'fields' defaults to -1, 'width' defaults to self.FIELD_BITS.
        """
        return smask_default(self.FIELD_BITS, value, fields, width)

    def gen_mmv_rotl(self, src, count, dest = None):
        """Rotate 'src' left by count 'bits', store result to 'dest'

        This generatates code for computing 
        'dest' = 'src' * 2**k  (mod p)
        for any fixed integer k. 'src' and 'dest' of type uint_mmv_t.
        They are interpreted as arrays of bit fields. The operation 
        is done on all bit fields simultaneously.
        """
        if dest is None:
            dest = src
        count = int(count) % self.P_BITS
        if count == 0:
            return
        s = """// Put %{dest} = %{src} rotated left by %{c} for all fields
%{dest} = (((%{src}) << %{c}) & %{smask: (-1 << {c}) & P})
       | (((%{src}) >> %{int:P_BITS-c}) & %{smask:range(c)});
"""
        return self.snippet(s, src = src, dest = dest, c=count)

    def gen_mmv_uint_spread(self, src, dest, value = None):
        """Spread bits of the integer 'src' to field array 'dest'

        'dest' is of type uint_mmv_t and interpreted an array of 
        bit fields. bit field i of 'dest' is set to -1 (mod p) if 
        bit i of 'src' is one and to 0 otherwise.

        Bit field i of 'dest' is set to  'value' if bit i of 'src' is 
        one. 0 <= 'value' <= self.p must hold. 'value' dafaults to
        self.p
        """ 
        msize = 1 << (self.INT_BITS >> 1)

        if value is None:
            value = self.P
        assert  0 <= value <= self.P
        if value == 0: 
            return "{dest} = 0;\n".format(dest = dest)
        s = """// Spread bits 0,...,{i} of {src} to the ({j}-bit long) fields
// of {dest}. A field of {dest} is set to {value} if its 
// corresponding bit in input {src} is one and to 0 otherwise.
""". format(i = self.INT_FIELDS - 1, value = hex(value), j = self.FIELD_BITS,
                      dest = dest, src= src)
        a, b = self.INT_FIELDS >> 1, self.INT_BITS
        while a:
            m0 = self.hex(self.smask(range(a), -1, b))
            m1 = self.hex(self.smask(range(a, 2*a), -1, b))
            b >>= 1
            s0 = "{dest} = ({src} & {m0}) + (({src} & {m1}) << {sh});\n"
            if len(m0 + m1) >= 17:
                s0 = "\n    + ".join(s0.split("+"))
            s += s0.format(dest=dest, src=src, m0=m0, m1=m1, sh=b-a)
            a >>= 1
            src = dest
        s += "%s *= %s;\n"  % (dest, value)
        s += "// Bit spreading done.\n"
        return s




class Mockup_MM_Op(MM_Op):
    def __init__(self, **kwds):
        self.p = p = int(kwds.get('p', 3))
        super(Mockup_MM_Op, self).__init__(p = p)
        old_tables = self.tables
        self.tables = {}
        self.tables.update(old_tables)
        self.tables["GENERATE_CODE"] = False



Tables = MM_Op
MockupTables = Mockup_MM_Op

