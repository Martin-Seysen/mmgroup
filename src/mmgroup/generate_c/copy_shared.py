import sys
import os
import subprocess
import re
import argparse
from shutil import copyfile





#############################################################################
# Copy shared libraries to build directory using BuildExtCmdObj object
#############################################################################


def has_extension(filename, extension_list):
    return any(filename.endswith(ext) for ext in extension_list)




def copy_shared_libs(build_lib = "", verbose = 1):
    """This is necessary for setuptools/distutils"""

    assert isinstance(build_lib, str)
    DIR = 'src/mmgroup'
    ABSDIR = os.path.abspath(DIR)
    if os.name in ["nt"]:
        extensions =  [".pyd", ".dll"]
    elif os.name in ["posix"]:
        extensions = [".so"]
    else:
        W = "don't know how do process shared libraries in a %s system"
        print("Warning:", W % os.name)
        extensions =  []

    files = os.listdir(ABSDIR)
    shared = [x for x in files if has_extension(x, extensions)]

    if build_lib.startswith('null') and not build_lib[4:5].isalnum():
        if verbose:
            print("No shared libraries to copy")
        return

    print("*** build_lib in setup.py =", build_lib)
    OUTDIR = os.path.abspath(os.path.join(build_lib, 'mmgroup'))

    for filename in shared:
        path = os.path.abspath(os.path.join(ABSDIR, filename))
        dest = os.path.join(OUTDIR, filename)
        if verbose:
            print("Copying %s to %s" % (path, dest))
        copyfile(path, dest)



if __name__ == "__main__":
    copy_shared_libs(sys.argv[1], verbose = 1)
