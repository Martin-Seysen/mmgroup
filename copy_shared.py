import sys


sys.path.append('src')

from mmgroup.generate_c.copy_shared import copy_shared_libs


if __name__ == "__main__":
    copy_shared_libs(sys.argv[1])


