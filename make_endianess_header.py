"""Detect endianess of local machine and write it into a header file"""

import sys
import os
import time
import argparse




def _comment_endianess(endianess, mockup = False):
    if 0 <= endianess <= 1:
        s = "Local machine is " + ['litte', 'Big'][endianess] 
        return s + " endian"
    else:
        explain = "has not been" if mockup else "could not be"
        return "Endianess %s detected for this machine!" % explain

def get_endianess(mockup = False, verbose = 0):
    if mockup:
        endianess = -1
    else:
        try:
            from mmgroup.mat24 import check_endianess
        except:
            time.sleep(0.1) # GitHub may be a bit slow here
            from mmgroup.mat24 import check_endianess
        endianess = check_endianess()
    if verbose:
        print(_comment_endianess(endianess, mockup))
    print("END", endianess, mockup)
    return endianess

s = r"""// This header has been generated automatically. Do not edit!
// It describes the endianess of the local machine.
// Do not copy this file to a different machine!
//
#ifndef {0}_H_INCLUDED
#define {0}_H_INCLUDED
// ENDIANESS is #defined to be 0 for little and 1 for Big endian.
// It is undefined if the endianess has not been detected.

"""

def write_header(h_file, mockup = False, verbose = 0):
    filename = os.path.splitext(os.path.split(h_file)[-1])[0]
    endianess = get_endianess(mockup)
    comment = _comment_endianess(endianess, mockup)
    f = open(h_file, "wt")
    f.write(s.format(filename.upper()) +  "// %s\n" % comment)
    if 0 <= endianess <= 1:
        f.write("#define ENDIANESS " + str(endianess) + "\n")
    f.write("#endif\n")
    f.close()
    if verbose:
        print(comment)


###########################################################################
# Command line interface
###########################################################################


def make_endianess_parser():
    description = ('Generate header file for the mmgroup project '
    'that defines the endianess of the local machine. '
    )

    # epilog = ("Some more documentation will follow here later."
    #)

    parser = argparse.ArgumentParser(
        prog = 'make_endianess_header',
        description=description,
        # epilog=epilog,
    )

    parser.add_argument('filename',  metavar='FILE', type=str,
         action = 'store', help='Name of generated header file')

    parser.add_argument('--mockup',
        default = False,
        action = 'store_true',
        help = 'Mockup output header for Sphinx')

    parser.add_argument('--library-path', nargs = '*', action='extend',
        metavar = 'PATHS', default = [],
        help = 'Set list of PATHS for finding shared libraries')

    parser.add_argument('--no-library-path', action = 'store_true',
        help=argparse.SUPPRESS)

    parser.add_argument('-v', '--verbose', action = 'store_true',
        help = 'Verbose operation')

    return parser




if __name__ == "__main__":
    sys.path.append('src')
    from mmgroup.generate_c import parse_set_shared_libraries
    sys.path.pop()

    # Set paths to shared libraries as given by command line args
    env_changed = parse_set_shared_libraries(sys.argv[1:])
    if env_changed:
        # Do the job in a subprocess, since the environment has changed
        import subprocess
        args = [
            sys.executable, '-m', 'make_endianess_header', '--no-library-path'
        ]
        # Same stupid Linux-like LD_LIBRARY_PATH stuff as
        # in module generate_code.py
        subprocess.check_call(args + sys.argv[1:])
    else:
        parsed_args = make_endianess_parser().parse_args(sys.argv[1:])
        write_header(parsed_args.filename, parsed_args.mockup,
            parsed_args.verbose)


