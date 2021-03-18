"""This is the *mmgroup* package for computing in the monster group.

In mathematics the *monster group* is the largest sporadic 
finite simple group.  

For documentation see:

https://mmgroup.readthedocs.io/en/latest/

The most important classes in this module are

  =================  =============================================
  ``class MMGroup``  Models an instance of the monster group.
  ``class MMSpace``  Models an instance of the 196884-dimensional
                     representation of the monster group
                     (modulo some small odd numbers).
  =================  =============================================

"""

import warnings 

try:
    from mmgroup.structures.autpl import  AutPL
    from mmgroup.structures.parity import Parity
    from mmgroup.structures.gcode import GCode, GcVector
    from mmgroup.structures.cocode import Cocode
    from mmgroup.structures.ploop import PLoopOne, PLoopOmega
    from mmgroup.structures.ploop import  PLoop, PLoopZ
    from mmgroup.structures.suboctad import Octad, SubOctad
    try:
        from mmgroup.mat24 import MAT24_ORDER
    except (ImportError, ModuleNotFoundError):
        # Import a pure python substitute for the mmgroup.mat24 extension
        # if the original extension has not been found. For background, 
        # see section 'Mock up C extension modules' in file
        # docs/source/conf.py.
        from mmgroup.dev.mat24.mat24_ref import Mat24
        MAT24_ORDER = Mat24.MAT24_ORDER
        del Mat24

except:
    w = "Extension mmgroup.mat24 not found, package not functional!"
    warnings.warn(w, UserWarning)



try:
    import mmgroup.mm
    import mmgroup.mm_group
    import mmgroup.mm_space
    from mmgroup.mm_group import MMGroup,  MMGroupWord, MM
    from mmgroup.mm import INT_BITS
except:
    w = "Extension mmgroup.mm not found, package not functional!"
    warnings.warn(w, UserWarning)


import mmgroup.generate_c


try:
    from mmgroup.mm_space import characteristics
    assert 3 in characteristics()
    from mmgroup.mm_space import MMSpace, MMSpaceVector, MMS
except:
    w = "Extension mmgroup.mm3 not found, package not functional!"
    warnings.warn(w, UserWarning)


try:
    from mmgroup.structures.xsp2_co1 import Xsp2_Co1
except:
    w = "Module 'mmgroup.structures.xsp2_co1' not found!"
    warnings.warn(w, UserWarning)
    
    