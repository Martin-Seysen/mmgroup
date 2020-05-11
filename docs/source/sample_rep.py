import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join("..", "..", "src")))

from mmgroup import MMGroup
from mmgroup import MMSpace
# Create an instance M of the monster group
M = MMGroup()
# Create representation space of M (modulo 3)
V = MMSpace(3, M)
# Create an element g of M 
g = M(('d', 0x123), ('p', 217821225))
# Create an vector v in V 
v = V(('T', 712, 13), (2, 'X', 0x345, 13))
# Let g operate on v by right multiplication
print(v * g)
