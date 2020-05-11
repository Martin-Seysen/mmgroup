

import sys
import os
import subprocess


if sys.platform.startswith("win32"):
    platform = "WIN32"
else:
    raise NotImplementedError("Platform %s not supported" % sys.platform)


def make_dll_win32(
    name, 
    sources, 
    libraries, 
    include_dirs,
    lib_path,
    dll_path,
    extra_args
    ):
    """Yet to be documented"""
    compile_args = [ "-c", "-Wall", "-DMS_WIN64"] + extra_args 
    for ipath in include_dirs:
        compile_args.append("-I" + os.path.realpath(ipath))
    link_args = ["-shared",  "-Wall", "-DMS_WIN64"] + extra_args 
    objects = []

    # Compile sources and add objects to list 'objects'
    for source in sources:
        args = compile_args[:]
        args.append(os.path.realpath(source))
        args.append("-o")
        obj = os.path.realpath(os.path.splitext(source)[0] + ".o")
        args.append(obj)
        objects.append(obj)
        print("gcc " + " ".join(args))
        subprocess.call("gcc " + " ".join(args)) 

    # link
    largs = link_args[:] + objects 
    for inc_dir in include_dirs:
        # same search path for include files an libraries
        largs.append("-L" + inc_dir)
    for library in libraries: 
        largs.append("-l" + library)
    dll_path =   os.path.join(dll_path, "%s.dll" % name)
    largs +=  ["-o",  dll_path ]
    def_path =   os.path.join(lib_path, "%s.def" % name)
    largs += [ "-Wl,--output-def," + def_path ]
    lib_path =   os.path.join(lib_path, "lib%s.lib" % name)
    largs += [ "-Wl,--out-implib," + lib_path ]
    print("gcc " + " ".join(largs))
    subprocess.call("gcc " + " ".join(largs)) 

    print("removing object files...")
    for obj in objects:
       os.remove(obj)
    print(dll_path + "\nhas been generated.\n")
    


def make_dll(
    name, 
    sources, 
    libraries = [], 
    include_dirs = [],
    lib_path = ".",
    dll_path = ".",
    extra_args = ""
    ):
    if platform == "WIN32":
        return make_dll_win32(name, sources, libraries,  include_dirs,
            lib_path, dll_path, extra_args )