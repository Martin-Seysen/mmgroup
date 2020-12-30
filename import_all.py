import sys 
import os
import time
sys.path.insert(0, os.path.abspath("src"))
#print(sys.path)

print("Importing python dlls..")
t_start = time.time()
import mmgroup 
from mmgroup import mat24
from mmgroup import generators
from mmgroup import mm
print("Importing python for specific characteristics..")
from mmgroup.mm_space import characteristics
pp = characteristics()
t = time.time() - t_start
print("done after %.2f seconds." % t)
s = "Representations of the Monster modulo %s are supported." 
print(s % ", ".join(map(str,pp)))

