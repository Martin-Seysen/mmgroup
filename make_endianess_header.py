"""Detect endianess of local machine and write it into a header file"""

import sys
import os
import time
sys.path.append('src')

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

if __name__ == "__main__":
    write_header(sys.argv[1], "--mockup" in sys.argv, verbose = 1)