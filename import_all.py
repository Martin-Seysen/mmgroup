import sys
sys.path.append("src")

print("Importing python dlls..")
import mmgroup.mat24
import mmgroup.mat24_xi
import mmgroup.mm


from mmgroup.mm_space import characteristics
characteristics()
print("done")

