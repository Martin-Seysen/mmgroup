import sys
import os
import subprocess
import re
import argparse
from shutil import copyfile


class MyArgumentParser(argparse.ArgumentParser):
    def convert_arg_line_to_args(self, arg_line):
        return arg_line.split()


def generate_linuxpatch_parser():
    description = ("Patch posix shared libraries with patchelf, so that"
    " everything works if all shared libraries are in same directory."
    )

    epilog = ("No action is performed on a non-posix system."
    )

    parser = MyArgumentParser(
        description=description, epilog=epilog,
        fromfile_prefix_chars='+',
    )

    parser.add_argument('--dir', 
        metavar = 'DIR',
        action = 'store',
        help = "Patch posix shared libraries (with extension '.so') "
               "in directory DIR."
    )

    parser.add_argument('-v', '--verbose', action = 'store_true',
        help = 'Verbose operation')

    return parser




############################################################################
# Patch posix shared libraries
############################################################################


def patch_shared(dir, verbose = 1):
    if os.name in ["nt"]:
        return     
    elif os.name in ["posix"]:
        return patch_shared_posix(dir, verbose)
    else:  
        W = "don't know how do process shared libraries in a %s system"
        print("Warning:", W % os.name)     
        return 




POSIX_MSG = """
Recommendation:

Install the 'patchelf' utility with the shell command:

sudo apt-get install patchelf

and rebuild the package!

"""

def patch_shared_posix(dir, verbose = 1):
    absdir = os.path.abspath(dir)
    shared = [x for x in os.listdir(absdir) if x.endswith('.so')] 

    if len(shared) == 0:
        print("Warning: No shared libraries found in directory %s" %
               absdir)
    print("Patching shared libraries in posix...")
    for filename in shared:
        path = os.path.abspath(os.path.join(absdir, filename))
        if verbose:
            print("patching " + path)
        args = ['patchelf', '--set-rpath', '$ORIGIN', path]
        try:
            subprocess.check_call(args)
        except:
            print('Executing', ' '.join(args))
            print('failed')
            print(POSIX_MSG)
            raise


#############################################################################
# Copy shared libraries to build directory using BuildExtCmdObj object
#############################################################################


def has_extension(filename, extension_list):
    return any(filename.endswith(ext) for ext in extension_list)


def copy_shared_libs(build_ext_cmd = None, verbose = 1):
    """This is necessary for setuptools/distutils"""

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

    if (not build_ext_cmd) or build_ext_cmd.inplace:
        return

    build_lib = build_ext_cmd.build_lib
    print("*** build_lib in setup.py =", build_lib)
    OUTDIR = os.path.abspath(DIR.replace('src', build_lib))


    for filename in shared:
        path = os.path.abspath(os.path.join(ABSDIR, filename))
        dest = os.path.join(OUTDIR, filename)
        if verbose:
            print("Copying %s to %s" % (path, dest))
        copyfile(path, dest)





if __name__ == "__main__":
    parser = generate_linuxpatch_parser()
    cmdline_args = parser.parse_args(sys.argv[1:])      
    patch_shared(cmdline_args.dir, cmdline_args.verbose)
