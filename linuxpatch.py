import sys


sys.path.append('src')

from mmgroup.generate_c.linuxpatch import patch


if __name__ == "__main__":
    patch(sys.argv[1:])


