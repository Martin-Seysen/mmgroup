import sys
sys.path.append('src')
from mmgroup.generate_c.shared_headers import copy_shared_all

if __name__ == "__main__":
    copy_shared_all(sys.argv[1:])



