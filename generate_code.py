"""Executable program for class generate_code.CodeGenerator"""

import sys
import argparse

py_path_parser = argparse.ArgumentParser(add_help=False)
py_path_parser.add_argument('--py-path',
        nargs = '*', action='extend', default = []
)

if __name__ == "__main__":
    py_path_args = py_path_parser.parse_known_args(sys.argv[1:])
    py_path = py_path_args[0].py_path + ['src'] 
    for i, path in enumerate(py_path):
        sys.path.insert(i, path)

    from mmgroup.generate_c import generate_code_parser, CodeGenerator
    parser = generate_code_parser()
    cmdline_args = parser.parse_args(sys.argv[1:])      

    cg = CodeGenerator(cmdline_args)
    #cg.activate_py_path()
    cg.import_tables()
    cg.display_args()
    cg.generate()
    #cg.deactivate_py_path()
    


