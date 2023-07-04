"""Auxiliary functions for writing tables into C program"""

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import sys
import re
import os
import argparse

from mmgroup.generate_c.generate_code import MyArgumentParser

from mmgroup.generate_c.generate_code import find_file
from mmgroup.generate_c.generate_code import open_for_write
from mmgroup.generate_c.make_pxd import pxd_from_h



def generate_pxd_parser():
    """Create parser for reading arguments for code generation 

    """
    description = ("Generate .pdx and .pxi file from header file"
    "for the mmgroup project using certain tables in the source file."
    )

    epilog = ("Some more documentation will follow here later."
    )

    parser = MyArgumentParser(
        description=description, epilog=epilog,
        fromfile_prefix_chars='+',
    )

    parser.add_argument('--pxd-out', 
        metavar='PXD',
        action='store', default = None,
        help="Set name of generated .pxd file to PXD."
    )

    parser.add_argument('--h-in', 
        metavar='HEADER',
        action='store', default = None,
        help = "Set input HEADER for generating .pxd file."
    )

    parser.add_argument('--pxd-in', 
        metavar='PXD',
        action='store', default = None,
        help="Set (optional) input .pxd file PXD for generating .pxd file."
    )

    parser.add_argument('--h-name', 
        metavar='HEADER',
        action='store', default = None,
        help="Set name of header to be stored in .pxd file to HEADER. "
             "Default is name of input header file."
    )

    parser.add_argument('--nogil', 
        action='store_true', default = false,
        help = "Optional, declare exported functions as 'nogil' "
        "when set."
    )
   
    parser.add_argument('--source-path', nargs = '*', action='extend',
        metavar = 'PATHS', default = [],
        help = 'Set list of PATHS for finding source files')

    parser.add_argument('--out-dir', 
        metavar = 'DIR', default = None,
        help = 'Set directory DIR for output file') 
 
    parser.add_argument('-v', '--verbose', action = 'store_true',
        help = 'Verbose operation')

    return parser



def check_args_parsed(args):
    ERR = False
    if not getattr(args, "pxd_out", None):
        ERR = "No output .pxd file specified"
    if not getattr(args, "h_in", None):
        ERR = "No input header file specified"
    if ERR:
        raise ValueError(ERR)

class pxdGenerator:
    """Yet to be documented"""
    def __init__(self, args):
        if isinstance(args, str):
            parser = generate_pxd_parser()
            self.s = parser.parse_args(args.split())
        elif isinstance(args, list):
            parser = generate_pxd_parser()
            self.s = parser.parse_args(args)
        elif isinstance(args, argparse.Namespace):
            self.s = args
        check_args_parsed(self.s)

    def find_file(self, filename):
        return find_file(self.s.source_path, filename)

    def generate(self):
        s = self.s
        h_in = os.path.norm_path(s.h_in)
        h_name = getattr(s, "h_name", None)
        h_name = os.path.norm_path(h_name) if h_name else h_in
        h_in = self.find_file(h_in)
        out_dir = getattr(s, "out_dir", None)
        pxd_out = open_for_write(out_dir, s.h_out)
        pxd_in = getattr(s, "pxd_in", None)
        pxd_in = self.find_file(pxd_in) if pxd_in else None
        nogil = getattr(s, "no_gil", False)
        pxd_from_h(pxd_out, h_in, pxd_in,  h_name, nogil)

 

def example():
    SAMPLE_ARGS = """ 
    """
    parser = generate_pcd_parser()
    s = parser.parse_args(SAMPLE_ARGS.split())
    print(s)
    print("\n")
    pg = pxdGenerator(s)
    print("")
    parser.print_help()


if __name__ == "__main__":
    example()

