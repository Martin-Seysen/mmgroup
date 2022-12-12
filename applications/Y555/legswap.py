r"""This little program is yet to be doumented.

We consider the automorphism ``s`` of the Y_555 diagram that exchanges
b_2 with b_3, c_2 with c_3, d_2 with d_3, etc. We will show that the
image ``img_s`` of ``s`` in the BiMonster is an element of the diagonal 
of the the direct product of the Monster with itself. Here the two
factors are 2A involutions in the Monster.

"""

try:
    import mmgroup
except (ImportError, ModuleNotFoundError):
    # get the mmgroup package from its inplace location it not installed
    import os
    import sys
    sys.path.append(os.path.join('..', '..', 'src'))
    import mmgroup


from mmgroup.bimm import BiMM, AutP3, AutP3_BiMM, P3_node
from mmgroup import MM

# Automorphism ``s`` is uniquely defiend by imagess of a, c1, c2, and c3
s = AutP3('a:a, c1:c1, c2:c3, c3:c2')
print("""Let 's' be the automorphism of the Y_555 graph that exchanges
b_2 with b_3, c_2 with c_3, d_2 with d_3, etc., and fixes a, b1, c1, etc.
Let 'img_s' be the image of 's' in the BiMonster.
""")
# Let ``img_s`` be the image of ``s`` in he BiMonster
img_s = AutP3_BiMM(s)
print("'img_s' is ", img_s)
# Decompose the image ``img_s``
m1, m2, e = img_s.decompose()
# Check that 'img_s' is in the direct square of the Monster 
assert e == 0
# Check that ``img_s`` is in the diagonal 
assert m1 == m2
print("'img_s' is in the diagonal of the direct square of the Monster.")
I1, _ =  m1.conjugate_involution()
assert I1 in [1,2]
class_ = '2A' if I1 == 1 else '2B'
msg = "The factors of 'img_s' are %s involutions in the Monster."
print(msg % class_) 
