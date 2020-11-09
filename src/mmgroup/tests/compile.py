import subprocess
import os


def compile_gcc(c_names, exe_name):
    args = ["gcc" ] + c_names + ["-o", exe_name]
    subprocess.check_output(args)

def compile_msvc(c_names, exe_name):
    args = ["cl" ] + c_names + ["/Fe" + exe_name]
    subprocess.check_output(args)

def compile_testprogramm(c_names, exe_name):
    """Compile a list of C programs to a single executable program
 
    ``c_names`` is a list of file names of C programs.
    ``exe_name`` is the file name of the executable program.
    """
    if os.name in ["nt"]:
        try:
            compile_msvc(c_names, exe_name)
            return
        except FileNotFoundError:
            pass
    compile_gcc(c_names, exe_name)
