r"""We deal with octads and certain subsets of octads called suboctads.

An *octad* is a  Golay code word of weight ``8``. The Golay code 
has ``759`` octads. A signed octad is a Parker loop element 
corresponding to an octad or a complement of an octad. Signed 
octads are used for indexing some unit vectors of the 
representation :math:`\rho` of the monster :math:`\mathbb{M}`.

Function ``Octad()`` returns a signed octad as an instance of 
class |PLoop|. Arguments of function ``Octad()`` are interpreted
in the same way as arguments of the constructor for class |PLoop|,
but the function raises `ValueError` if the result is not a
(possibly complemented) signed octad.

We use some lexical order for numbering the ``759`` octads,
which we do not describe in detail.  
Due to the well-known properties of the Golay code we can create 
a random octad and display its number as follows:

.. code-block:: python
    
        from random import sample
        from mmgroup import Octad
        # Select 5 random integers between 0 and 23
        int_list = sample(range(24), 5)
        # Complete these 5 integers to a (unique) octad o
        o = Octad(int_list)
        # Display octad o and its number
        print("Octad", o.bit_list, "has number", o.octad)


A *suboctad* of an octad ``o`` is an element of the Golay cocode
:math:`\mathcal{C}^*` of even parity which can be represented
as a subset of ``o``. Each octad has ``64`` suboctads.
A pair (octad, suboctad) is used for indexing some basis
vectors in the representation :math:`\rho`. function ``SubOctad``
creates a suboctad as an instance of class |XLeech2| from such a 
pair. Function ``SubOctad`` takes two arguments ``octad`` and 
``suboctad``. The first argument evaluates to a signed octad as 
in function ``Octad()``. The second argument evaluates to a 
suboctad, see  function ``SubOctad`` for details.

The raison d'etre of a  suboctad is indexing a basis vector in
the representation  :math:`\rho`. For this purpose we need a pair 
of integers refering to the octad and the suboctad. For an instance 
``so`` of class |XLeech2| that pair is given as the pair of the
last two integers in ``so.vector_tuple()``.

For an octad ``o`` and a suboctad ``s`` precisely one of the
pairs ``(o,s)`` and ``(~o,s)`` is actually used as an index for 
:math:`\rho`. The pair ``(~o,s)`` is used if the cocode element 
corresponding to ``s`` has minimum weight ``2``; the pair 
``(o,s)`` is used if that element has minimum weight ``0`` or 
``4``. See :cite:`Con85` or :cite:`Seysen20` for background. In 
this  implementation both, ``(o,s)`` and ``(~o,s)``, refer to the 
same basis vector of :math:`\rho`.
"""



from functools import reduce
from operator import __xor__
from numbers import Integral, Number
from random import randint
import warnings




from mmgroup.structures.abstract_group import AbstractGroupWord
from mmgroup.structures.parity import Parity
from mmgroup.structures.parse_atoms import ihex

from mmgroup.structures.autpl import AutPL, AutPlGroup

from mmgroup.structures.gcode import GCode, GcVector
from mmgroup.structures.cocode import Cocode
from mmgroup.structures.ploop import PLoop


ERR_RAND = "Illegal string for constricting type %s element" 

ERR_DIV4 = "%s object may be divided by 2 or 4 only"
ERR_DIV2 = "%s object may be divided by 2 only"



#######################################################################
# Import derived classed
#######################################################################

import_pending = True

def complete_import():
    """Internal function of this module

    If you need any of the objects declared above, do the following:

    if import_pending:
        complete_import()
    """
    global import_pending, mat24, XLeech2
    from mmgroup import mat24
    from mmgroup.structures.xleech2 import XLeech2
    import_pending = False



#######################################################################
# Auxiliary functions
#######################################################################



def as_vector24(data):
    """Convert list of bit positions to vector in ``GF(2)**24``"""
    try:
        vector = reduce(__xor__, (1 << x for x in data), 0)
    except:
        err = "List of integers >= 0 required for a 24-bit vector"
        raise TypeError(err)
    if vector >= 0x1000000:
        raise ValueError("List entry too large for a 24-bit vector")
    return vector


def str_vector24(v):
    """Return the 24-bit vector ``v`` as a string of binary digits

    Here ``v`` must be an integer with the coefficient of the
    ``i``-th bit coded in the bit with valence ``2**i``.
    Coefficients  ``i >= 24`` are ignored.
    """
    l = ["01"[(v >> i) & 1] for i in range(24)]
    l4 = ("".join(l[i:i+4]) for i in range(0, 24, 4))
    return " ".join(l4)

def str_basis(text, letter, basis):
   s = text + "\n"
   for i, v in enumerate(basis):
       s += "  %3s: %s\n" % (letter + str(i), str_vector24(v)) 
   return s


#######################################################################
# Function Octad
#######################################################################



def Octad(octad):
    r"""Return a (signed) octad as an element of the Parker loop. 

    :param octad:  the value of the octad to be returned. 
    :type octad: see table below for legal types

    :return:  a (signed) octad
    :rtype:   an instance of class |PLoop|

    :raise:
        * TypeError if ``type(octad)`` is not in the table given above.
        * ValueError if ``octad`` cannot be converted to a
          (possibly negated and complemented) octad.


    Depending on its type parameter **octad** is  interpreted as follows:

    .. table:: Legal types for parameter ``octad``
      :widths: 20 80

      ===================== ================================================
      type                  Evaluates to
      ===================== ================================================
      ``int``               Here the (positive) octad with  number 
                            ``octad`` is returned. There are 759 octads
                            in the Golay code. 
                            So ``0 <= octad < 759`` must hold.
                            

      ``list`` of ``int``   Such a list is converted to a Golay code word,
                            see class |GCode|, and the corresponding 
                            (positive) Parker loop element is returned.

      class |GCode|         The corresponding 
                            (positive) Parker loop element is returned. 

      class |PLoop|         A deep copy of parameter ``octad`` is returned.

      class |GcVector|      This is converted to a Golay code word,
                            see class |GCode|, and the corresponding 
                            (positive) Parker loop element is returned.

      class |XLeech2|       The *octad* part of the vector ``octad``  
                            is  returned. 

      ``str``               Create random element depending on the string
                             | ``'r'``: Create arbitrary octad

      ===================== ================================================

    A complement of an octad is also accepted; then the corresponding 
    Parker loop element is retured. The function raises ValueError
    if  parameter ``octad`` does not evaluate to an octad or a
    complement of an octad.      
    """
    if import_pending:
        complete_import()
    if isinstance(octad, Integral):
        if  not 0 <= octad < 759:
            raise ValueError("Bad octad number")
        return PLoop(mat24.octad_to_gcode(octad) & 0xfff)
    elif isinstance(octad, str):
        return PLoop(mat24.octad_to_gcode(randint(0, 758)))
    elif isinstance(octad, PLoop):
        _ = octad.octad
        return octad
    elif isinstance(octad, GCode):
        _ = octad.octad
        return PLoop(octad.ord)
    else:
        v =  PLoop(octad)
        _ = v.octad
        return v


#######################################################################
# Function ectad_entries
#######################################################################

def octad_entries(octad):
    r"""Return list of entries of an octad

    Here parameter ``octad`` describes an octad as in function
    ``Octad``. That octad is a sum of eight basis vectors of
    :math:`\mbox{GF}_2^{24}`, with basis vectors numbered from 0 to 23.
    The function returns the list ``[b0,...,b7]`` of the indices of
    these eight basis vectors in a certain order that may differ from
    the natural order.

    We use these indices for numbering the suboctads of that octad
    as described in the documentation of function ``SubOctad``.

    .. Caution::
      The order of the basis vectors returned by this function
      may be changed in future versions!
    """
    if import_pending:
        complete_import()
    o = Octad(octad).octad
    return mat24.octad_entries(o)

#######################################################################
# Function  SubOctad
#######################################################################


def SubOctad(octad, suboctad = 0):
    r"""Creates a suboctad as instance of class |XLeech2|

    :param octad:

      The first component of the pair *(octad, suboctad)* to be 
      created.

    :type octad: same as in function
      :py:class:`~mmgroup.Octad`

    :param suboctad:
    
      The second component of the pair 
      *(octad, suboctad)* to be created.

    :type suboctad: see table below for legal types


    :return: 

     The suboctad given by the pair *(octad, suboctad)* 

    :rtype: an instance of class |XLeech2|

    :raise:
        * TypeError if one of the arguments *octad* or *suboctad* 
          is not of correct type.
        * ValueError if argument *octad* or *suboctad* does not
          evaluate to an octad or to a correct suboctad,
          respectively.
    

    A *suboctad* is an even cocode element that can be represented
    as a subset of the octad given by the argument *octad*. 

    The raison d'etre of function ``SubOctad`` is that pairs
    *(octad, suboctad)* are used for indexing vectors in the
    representation of the monster group. Here we want to 
    number the octads from ``0`` to ``758`` and the suboctads
    form ``0`` to ``63``, depending on the octad. Note that
    every octad has ``64``  suboctads.

    Depending on its type parameter **suboctad** is  interpreted as follows:

    .. table:: Legal types for parameter ``suboctad``
      :widths: 20 80

      ===================== ================================================
      type                  Evaluates to
      ===================== ================================================
      ``int``               Here the suboctad with the number given 
                            in the argument is 
                            taken.  That numbering depends on the octad 
                            given in   the argument ``octad``. 
                            ``0 <= suboctad < 64`` must hold.                           

      ``list`` of ``int``   Such a list is converted to a bit vector
                            as in class |GcVector|,
                            and the cocode element corresponding to that
                            bit vector is taken.

       class |GCode|        The intersection of the octad given as the 
                            first argument and the Golay code word given
                            as the second argument is taken. 
  
       class |GcVector|     This is converted to a cocode element,
                            see class |Cocode|, and that cocode element 
                            is taken.

       class |Cocode|       That cocode element is taken as the suboctad.


      ``str``               Create random element depending on the string
                             | ``'r'``: Create arbitrary suboctad
      ===================== ================================================

    The numbering of the suboctads

    Suboctads are numbered for 0 to 63. Let ``[b0, b1,..., b7]`` be the 
    bits set in the octad of the pair ``(octad, suboctad)``. Here we
    assume the indices of the basis vectors making up the octad are
    ordered as returned by function ``octad_entries``, when applied
    to the octad.

    The following table shows the suboctad numbers for some suboctads
    given as cocode elements. More suboctad numbers can be obtained by
    combining suboctads and their corresponding numbers  with ``XOR``.

    .. table:: Suboctad numbers of some cocode elements
      :widths: 16 16 16 16 18 18

      =========== =========== =========== =========== =========== =========== 
      ``[b0,b1]`` ``[b0,b2]`` ``[b0,b3]`` ``[b0,b4]`` ``[b0,b5]`` ``[b0,b6]``
      ``s  = 1``  ``s  = 2``  ``s  = 4``  ``s  = 8``  ``s = 16``  ``s = 32``
      =========== =========== =========== =========== =========== =========== 

    E.g. ``[b0, b5, b6, b7]`` is equivalent to ``[b1, b2, b3, b4]`` 
    modulo the Golay code and has number ``s = 1 ^ 2 ^ 4 ^ 8 = 15``.

    """
    if import_pending:
        complete_import()
    ploop = Octad(octad)
    gcode = ploop.value & 0xfff
    if isinstance(suboctad, str):
        suboctad_ = randint(0, 63)
    elif isinstance(suboctad, Integral):
        suboctad_ = suboctad & 0x3f
    elif isinstance(suboctad, GCode):
        value = mat24.ploop_cap(gcode, suboctad.value)
        suboctad_ = mat24.cocode_to_suboctad(value, gcode) & 0x3f
    else:
        value = Cocode(suboctad).cocode 
        suboctad_ = mat24.cocode_to_suboctad(value, gcode) & 0x3f
    cocode = mat24.suboctad_to_cocode(suboctad_, ploop.octad)
    result = XLeech2(ploop, cocode)
    subtype =  result.xsubtype
    assert subtype in [0x22, 0x42], (hex(subtype), hex(ploop), hex(cocode), hex(result.value))
    # if necessary, complement octad to make Leech vector short
    if subtype == 0x42:
        result.value  ^= 0x800000
    return result




