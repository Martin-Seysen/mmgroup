import sys
import os
import argparse
import shutil




def build_parser():
    description = ('Copy headers for shared libraries into wheel'
    )

    # epilog = ("Some more documentation will follow here later."
    #)

    parser = argparse.ArgumentParser(
        prog = 'shared_shared',
        description=description,  
        # epilog=epilog,
    )


    parser.add_argument('--header-path', 
        nargs = '*',  metavar='PATHS',
        action = 'extend', default = [], 
        help = 'Set list of PATHS for finding .h files to be copied.')

    parser.add_argument('--header-dir', 
        metavar='DIR', action='store', default = "",
        help = 'Set directory DIR for storing .h files.'
    )

    parser.add_argument('--lib-path', 
        nargs = '*',  metavar='PATHS',
        action = 'extend', default = [], 
        help = 'Set list of PATHS for finding import libraries to be copied. '
               'Used in Windows, ignored in Unix-like systems' 
        )

    parser.add_argument('--lib-dir', 
        metavar='DIR',
        action='store', default = "",
        help = 'Set directory DIR for storing import  libraries. '
               'Used in Windows, ignored in Unix-like systems' 
    )

    return parser



def copy_directory(source_dir, extension, dest_dir):
    os.makedirs(dest_dir, exist_ok=True)
    files = [f.name for f in os.scandir(source_dir) if f.is_file()]
    for fname in files:
        if os.path.splitext(fname)[1] == extension:
            src_path = os.path.join(source_dir, fname)
            dest_path = os.path.join(dest_dir, fname)
            shutil.copyfile(src_path, dest_path)
    
def copy_headers(cmdline_args):
    if cmdline_args.header_dir:
        dest_dir = cmdline_args.header_dir 
        for source_dir in cmdline_args.header_path:
            copy_directory(source_dir, '.h', dest_dir) 

def copy_libs(cmdline_args):
    if os.name == 'nt' and cmdline_args.lib_dir:
        dest_dir = cmdline_args.lib_dir 
        for source_dir in cmdline_args.lib_path:
            copy_directory(source_dir, '.lib', dest_dir) 


def copy_shared_all(args):
    parser = build_parser()
    cmdline_args = parser.parse_args(args)
    copy_headers(cmdline_args)
    copy_libs(cmdline_args)


if __name__ == "__main__":
    copy_shared_all(sys.argv[1:])


  
