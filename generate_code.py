"""Executable program for class generate_code.CodeGenerator"""

import sys
import argparse

path_parser = argparse.ArgumentParser(add_help=False)
path_parser.add_argument('--library-path',
        nargs = '*', action='extend',  default = []
)
path_parser.add_argument('--no-library-path', action = 'store_true')


if __name__ == "__main__":
    sys.path.append('src')
    from mmgroup.generate_c import generate_code_parser, CodeGenerator
    from mmgroup.generate_c.generate_code import set_shared_libraries
    sys.path.pop()

    path_args = path_parser.parse_known_args(sys.argv[1:])
    env_changed = set_shared_libraries(path_args[0])
    if env_changed:
        # Do the job in a subprocess, since the environment has changed
        import subprocess
        args = [
            sys.executable, '-m', 'generate_code', '--no-library-path'
        ]
        subprocess.check_call(args + sys.argv[1:])
    else:
        parser = generate_code_parser()
        cmdline_args = parser.parse_args(sys.argv[1:])

        cg = CodeGenerator(cmdline_args)
        cg.activate_py_path()
        cg.import_tables()
        cg.display_args()
        cg.generate()
        #cg.deactivate_py_path()
