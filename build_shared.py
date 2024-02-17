import sys

sys.path.append('src')

from mmgroup.generate_c.build_shared import build


if __name__ == "__main__":
    build(sys.argv[1:])


  
