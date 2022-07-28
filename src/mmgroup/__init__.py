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


try:
    # Try importing the fast C function
    from mmgroup.mat24 import MAT24_ORDER 
except (ImportError, ModuleNotFoundError):
    # Import a pure python substitute for the mmgroup.mat24 extension
    # if the original extension has not been found. For background, 
    # see section 'Mock up C extension modules' in file
    # docs/source/conf.py.
    from mmgroup.dev.mat24.mat24_ref import  Mat24
    mat24 = Mat24    
    MAT24_ORDER = Mat24.MAT24_ORDER
    del Mat24
    w = "Extension mmgroup.mat24 not found, package not functional!"
    warnings.warn(w, UserWarning)


try:
    from mmgroup.structures.parity import Parity
    from mmgroup.structures.gcode import GCode, GcVector
    from mmgroup.structures.cocode import Cocode
    from mmgroup.structures.ploop import PLoopOne, PLoopOmega
    from mmgroup.structures.ploop import  PLoop, PLoopZ
    from mmgroup.structures.suboctad import Octad, SubOctad
    from mmgroup.structures.autpl import  AutPL
    from mmgroup.structures.xleech2 import XLeech2
except:
    w = "Extension mmgroup.mat24 not found, package not functional!"
    warnings.warn(w, UserWarning)


try:
    import mmgroup.mm
    from mmgroup.mm import INT_BITS
except:
    w = "Extension mmgroup.mm not found, package not functional!"
    warnings.warn(w, UserWarning)
    

try:
    from mmgroup.structures.mm0_group import MM0Group,  MM0
except:
    w = "Class mmgroup.structures.MM0 not found, package not functional!"
    warnings.warn(w, UserWarning)



try:
    import mmgroup.mm_group
    from mmgroup.mm_group import MMGroup, MM
except:
    w = "Class mmgroup.MM not found, package not functional!"
    warnings.warn(w, UserWarning)



import mmgroup.generate_c


try:
    import mmgroup.mm_space
    from mmgroup.mm_space import characteristics
    assert 3 in characteristics()
    from mmgroup.mm_space import MMSpace, MMVector, MMV, order_vector
except:
    w = "Extension mmgroup.mm3 not found, package not functional!"
    warnings.warn(w, UserWarning)

try:
    from mmgroup.structures.xsp2_co1 import Xsp2_Co1
except:
    w = "Module 'mmgroup.structures.xsp2_co1' not found!"
    warnings.warn(w, UserWarning)
    
    
try:
    from mmgroup.mm_crt_space import  MMVectorCRT
except:
    pass


