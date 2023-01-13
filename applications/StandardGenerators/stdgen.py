r"""This application computes standard generators of the Monster.

The standard generators ``(a, b)`` of the Monster satisfy the following
properties:

``a`` is in class 2A, ``b`` is in class 3B and ``ab`` has order 29, 

see

https://brauer.maths.qmul.ac.uk/Atlas/v3/spor/M/

Function ``standard_generators`` in file ``stdgen.py`` in this 
application computes standard generators of the Monster (if this has not 
yet been done) and stores them in the file ``standard_generators.py``. 
"""

try:
    import mmgroup
except (ImportError, ModuleNotFoundError):
    # get the mmgroup package from its inplace location if not installed
    import os
    import sys
    sys.path.append(os.path.join('..', '..', 'src'))
    import mmgroup



from mmgroup import MM, Xsp2_Co1

def find_3B_element():
    """Return random element of class 3B of the Monster.

    The function return an element of class 3B of the Monster
    as an instance of class MM.

    The character of an element of class 3B in the rep of the
    Monster of dimension 196883 is 53.
    """
    # Compute in subgroup ``G_x0`` since this is faster
    while(1):
        b0 = Xsp2_Co1('r','G_x0') # random element of G_x0
        o = b0.order()            # order of that element
        if o % 3 == 0:            # Is order divisible by 3?
            b = b0 ** (o//3)      # Power up to an element b of order 3
            if b.chi_G_x0()[0] == 53:
                return MM(b)      # return b if it has character 53

def make_standard_generators(verbose = 1):
    """Find standard generators (a,b) of the Monster.

    The function finds standard generators (a,b) of the Monster and
    returns them as a pair of instances of class MM.
    """
    n_trials = 0
    a = MM('d', [2,3])              # This is generator a.
    assert a.order() == 2           # Check that a has order 2
    assert a.chi_G_x0()[0] == 4371  # and character 4371,
                                    # as required for a 2A element.
    b_conj = find_3B_element()      # b_conj is a 3B element.
    while 1:
        n_trials += 1
        if verbose:
            print("trial %d   \r" %  n_trials, end = "")
        b = b_conj ** MM('r', 'M')  # b is a random 3B element.
        if (a * b).order(29) == 29: # Done if  a*b  has order 29.
            if verbose:
                 s = "\n%s trials required for finding standard generators"
            print(s % n_trials)
            return a, b
       
   
def store_standard_generators(a, b):
    """Store standard generators of the Monster in a python script.

    The function stores the generators (a ,b), which should be 
    instances of class MM as computed by function
    ``make_standard_generators``, in the python script
    ``standard_generators.py``.
    """
    FILENAME = "standard_generators.py"
    output = """# This file has been computed automatically, do not change!
# Standard generators of the Monster group.
from stdgen import MM
a = MM(\"%s")
b = MM(\"%s\")
""" % (a, b)
    f = open(FILENAME, "wt")
    print(output, file = f)
    f.close()

     
def standard_generators():
    """Returns standard generators (a,b) of the Monster.

    The function returns standard generators (a,b) of the Monster
    as a pair of instances of class MM. For the definition of 
    standard generators see

    https://brauer.maths.qmul.ac.uk/Atlas/v3/spor/M/

    Warning:
    Standard generators are version dependent, since their 
    computation is based on a random generator!

    The function returns the precomputed generators in file 
    ``standard_generators`` if present. Otherwise it computes
    standard generators and stores them in that file.    
    """
    try:
        from standard_generators import a, b
        return a, b
    except (ImportError, ModuleNotFoundError):
        a, b = make_standard_generators()
        store_standard_generators(a, b)
        return a, b


if __name__ == '__main__':
     a, b = standard_generators()
     print("Standard generators a, b of the Monster group")
     print("a =", a)
     print("b =", b)
     
