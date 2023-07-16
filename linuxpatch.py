import sys
import os
import subprocess
import re
import argparse


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


if __name__ == "__main__":
    parser = generate_linuxpatch_parser()
    cmdline_args = parser.parse_args(sys.argv[1:])      
    patch_shared(cmdline_args.dir, cmdline_args.verbose)
