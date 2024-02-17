"""Executable program for class generate_code.CodeGenerator"""

import sys


if __name__ == "__main__":
    sys.path.append('src')
    from mmgroup.generate_c import CodeGenerator
    from mmgroup.generate_c import parse_set_shared_libraries
    sys.path.pop()

    # Set paths to shared libraries as given by command line args
    env_changed = parse_set_shared_libraries(sys.argv[1:])
    if env_changed:
        # Do the job in a subprocess, since the environment has changed
        import subprocess
        args = [
            sys.executable, '-m', 'generate_code', '--no-library-path'
        ]
        subprocess.check_call(args + sys.argv[1:])
    else:
        cg = CodeGenerator(sys.argv[1:])
        cg.activate_py_path()
        cg.import_tables()
        cg.display_args()
        cg.generate()
        #cg.deactivate_py_path()
