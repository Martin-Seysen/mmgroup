"""This is the *mmgroup* package for computing in the monster group.

In mathematics the *monster group* is the largest sporadic 
finite simple group.  

The most important classes in module ``mmgroup`` are

  ===================  ===========================================
  ``class MM``         Models an element of the monster group.
  ``class MMVector``   Models a vector in a 196884-dimensional
                       representation of the monster group
                       (modulo a small odd number).
  ===================  ===========================================

For documentation see:

https://mmgroup.readthedocs.io/en/latest/
"""

import warnings 


_user_has_been_warned = False
def _warn(message):
    global _user_has_been_warned
    if not _user_has_been_warned:
        warnings.warn(message, UserWarning, stacklevel=2)
        _user_has_been_warned = True


try:
    # Try importing the fast C function
    from mmgroup import mat24
    from mmgroup.mat24 import MAT24_ORDER 
    assert type(MAT24_ORDER) == int
    GCODE_BASIS = mat24.basis[12:24]
    COCODE_BASIS = mat24.basis[:12]
except (ImportError, ModuleNotFoundError, AssertionError):
    # Import a pure python substitute for the mmgroup.mat24 extension
    # if the original extension has not been found. For background, 
    # see section 'Mock up C extension modules' in file
    # docs/source/conf.py.
    from mmgroup.dev.mat24.mat24_ref import  Mat24
    mat24 = Mat24    
    MAT24_ORDER = Mat24.MAT24_ORDER
    GCODE_BASIS = mat24.basis[12:24]
    COCODE_BASIS = mat24.basis[:12]
    del Mat24
    w = "Extension mmgroup.mat24 not found, package not functional!"
    _warn(w)


try:
    _m = "parity"
    from mmgroup.structures.parity import Parity
    _m = "gcode"
    from mmgroup.structures.gcode import GCode, GcVector
    _m = "cocode"
    from mmgroup.structures.cocode import Cocode
    _m = "ploop"
    from mmgroup.structures.ploop import PLoopOne, PLoopOmega
    from mmgroup.structures.ploop import  PLoop, PLoopZ
    _m = "suboctad"
    from mmgroup.structures.suboctad import Octad, SubOctad
    from mmgroup.structures.suboctad import octad_entries
    _m = "autpl"
    from mmgroup.structures.autpl import  AutPL
    _m = "xleech2"
    from mmgroup.structures.xleech2 import XLeech2
    from mmgroup.structures.xleech2 import leech2_orbits_raw
except:
    #import traceback
    #traceback.print_stack()
    w = "The %s package is not functional!"
    _warn(w % _m)
del _m



try:
    import mmgroup.mm_op
    from mmgroup.mm_op import INT_BITS
except:
    w = "Extension mmgroup.mm_op not found, package not functional!"
    _warn(w)
    

try:
    from mmgroup.structures.mm0_group import MM0Group,  MM0
except:
    w = "Class mmgroup.structures.MM0 not found, package not functional!"
    _warn(w)



try:
    import mmgroup.mm_group
    from mmgroup.mm_group import MMGroup, MM, MM_from_int
except:
    w = "Class mmgroup.MM not found, package not functional!"
    _warn(w)



try:
    import mmgroup.mm_space
    from mmgroup.mm_space import characteristics
    #assert 3 in characteristics()
    from mmgroup.mm_space import MMSpace, MMVector, MMV, order_vector
    from mmgroup.mm_space import mmv_scalprod
except:
    w = "Extension mmgroup.mm_op not found, package not functional!"
    _warn(w)

try:
    from mmgroup.structures.xsp2_co1 import Xsp2_Co1
except:
    w = "Module 'mmgroup.structures.xsp2_co1' not found!"
    _warn(w)
    
    
try:
    from mmgroup.mm_crt_space import  MMVectorCRT
except:
    w = "Module 'mmgroup.mm_crt_space' not found!"
    _warn(w)


try:
    from mmgroup.generators import gen_rng_seed_init
    gen_rng_seed_init()
except:
    w =  "The random generator of this package is not functional!"
    _warn(w)

