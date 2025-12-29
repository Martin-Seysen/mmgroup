import sys
import re
import os
import argparse
import subprocess
import platform


from mmgroup.generate_c.parallel_processes import SimpleProcessWorkers




def shared_lib_name(name, mode, static=False, os_name=None, pymod=False, flat=False):
    """Adjust name of a (shared) library

    Here parameter ``name`` is a name of a (shared) library to be
    adjusted to various situations in the build process. This name
    may be prefixed by some path information, but the name itself
    should not contain any suffixes, e.g. '.a', '.so', '.lib', '.dll'.
    The part of parameter ``name`` indicating the library name
    should not contain any prefix strings, e.g. 'lib'.

    The function adjusts the name of the library as configured by the
    other parameters of the function and returns adjusted name as a
    string. If ``name`` is a list of strings then the function returns
    the list of adjusted strings.

    Parameter ``mode`` describes the part of the build process
    where the adjusted name is to be used:

    mode = 'load':      The file name name of the (shared) library
                        including all prefexes and extensions.

    mode = 'build':     The name of the library to be linked to
                        a program using that library. In Unix-like
                        systems this is the (shared) library itsalf.
                        In case of a Windows DLL this is an import
                        library generated together with the DLL.

    mode = 'build_ext'  Similar to mode 'build'; to be used in
                        ``setuptools`` instead of mode 'build'.

    Parameter ``static`` should be True in case of a static library
    and False (default) in case of a shared library.

    Parameter ``os_name`` defaults to the string returned by function
    ``os.name`` in the ``os`` package. It may be set to model the
    behaviour of a different os.

    If parameter ``pymod`` is True then ``name`` is parsed as a python
    module separated by '.' characters. Default is False.

    If parameter ``flat`` is True then any information in parameter
    ``name`` refering to a directory is removed. Default is False.
    """
    if not os_name:
        os_name = os.name
    if isinstance(name, list):
        return [shared_lib_name(x, mode, static, os_name, pymod, dir)
                 for x in name]
    if pymod:
        name = re.sub('\.', '/', name)
    path = os.path.split(name)
    if len(path) == 0:
        raise ValueError("No library name specified")
    path, name = list(path[:-1]), path[-1]
    if flat:
        path = []
    platform_found = new_name = False
    if os.name == "posix":
        platform_found = True
        if mode in ["build", "build_ext"]:
            new_name = name
        elif mode in ["load"]:
            suffix = ".a" if static else ".so"
            new_name =  "lib" + name + suffix
    elif os.name == "nt":
        platform_found = True
        if static:
            if mode == "build_ext":
                new_name = name
            elif mode in ["build", "load"]:
                new_name = name + ".lib"
        else:
           if mode == "build_ext":
               new_name = "lib" + name
           if mode == "build":
               new_name =  "lib" + name + ".lib"
           if mode == "load":
               new_name =  name + ".dll"

    if not platform_found:
        err = "Don't know the suffix for a %slibrary on platform %s"
        sh = '' if static else 'shared '
        raise ValueError(err % (sh, os_name))
    if not new_name:
        err = "Illegal mode %s in function shared_lib_name"
        raise ValueError(err % mode)
    new_name = os.path.normpath(os.path.join(*(path + [new_name])))
    #print('%s' % new_name)
    return new_name



def default_compiler():
    if os.name == "posix":
        if platform.uname().system == "Darwin":
            return 'darwin'
        return  'unix'
    elif os.name == "nt":
        return 'msvc'
    else:
        return 'unknown'






def build_shared_lib_parser():
    description = ('Generate a shared library from a collection of C sources'
    )

    # epilog = ("Some more documentation will follow here later."
    #)

    parser = argparse.ArgumentParser(
        prog = 'build_shared',
        description=description,  
        # epilog=epilog,
    )

    parser.add_argument('--name', required=True,
        action='store', default = None, metavar='MAME',
        help = "Set NAME of library to be generated. "
        "NAME should have no extension and no prefix such as 'lib'."
    )

    parser.add_argument('--sources', 
        nargs = '*',  metavar='SOURCES',
        action = 'extend', default = [], 
        help = "Add list SOURCES of to list of source files. "
        "Each source file should be a .c file with extension '.c'."
    )

    parser.add_argument('--source-dir', 
        metavar='DIR',
        action='store', default = "",
        help="Set root directory DIR for all source files."
    )

    parser.add_argument('--include-path',
        nargs = '*',  metavar='PATHS',
        action = 'extend', default = [], 
        help = 'Set list of PATHS for finding .h files to be used.')

    parser.add_argument('--library-path',
        nargs = '*',  metavar='PATHS',
        action = 'extend', default = [], 
        help = 'Set list of PATHS for finding libraries to be used.')

    parser.add_argument('--define',
        nargs = '*',  metavar='VARS',
        action = 'extend', default = [], 
        help = 'Define variables given by the list VARS for the compiler. '
        'Each variable must be a string VAR or VAR=VALUE.')

    parser.add_argument('--undef',
        nargs = '*',  metavar='VARS',
        action = 'extend', default = [], 
        help = 'Undefine variables given by the list VARS for the compiler.'
        )

    parser.add_argument('--libraries',
        nargs = '*',  metavar='LIBS',
        action = 'extend', default = [], 
        help = "Search libraries with names in the list given by LIBS. "
               "This corresponds to the gcc option '-l'."
        )

    parser.add_argument('--library-dir', 
        metavar='DIR', action='store', default = "",
        help = 'Set directory DIR for storing object files and static libraries.'
    )

    parser.add_argument('--shared-dir', 
        metavar='DIR',
        action='store', default = "",
        help = 'Set directory DIR for storing shared libraries.'
    )

    parser.add_argument('--lflags',
        metavar='FLAGS', action = 'store', default = '',
        help = "Add extra arguments FLAGS for the linker. E.g. "
               "'--lflags=-g' adds arguments '-g'."
        )

    parser.add_argument('--static', nargs='?', const=1, default=0, type=int,
         help = 'Create static instead of shared library. '
           'Optional argument STATIC may be 0 (=shared) or 1 (=static).'
        )

    parser.add_argument('--n',
        type=int,  default = 1, metavar = 'N',
        help = 'Use N parallel processes for compiling.'
        )
 
    parser.add_argument('--compiler',
        action = 'store', default = default_compiler(), metavar = 'C',
        help = "Specify name of default compiler. "
               "C must be 'unix', 'darwin', 'msvc', or 'mingw32'."
        )

    parser.add_argument('--cflags',
        metavar='FLAGS', action = 'store', default = '',
        help = "Add extra arguments FLAGS for the compiler. E.g. "
               "'--cflags=-c,foo=bar' adds arguments '-c' and 'foo=bar'."
        )

    parser.add_argument('--rpath',
        nargs = '*',  metavar='PATH',
        action = 'extend', default = [],
        help = "Set path PATH[:PATH] in unix-like os shared library with linker "
               "option -rpath. Ingnored in a non unix-like os."
        )

    parser.add_argument('--mockup', nargs='?', const=1, default=0, type=int,
         help = "Do nothing when set; for compatibility with the corresponding "
           "option in the 'generate_code.py' tool."
        )

    parser.add_argument('--display',
        action = 'store',  metavar = 'MODE',  default = None,
        help = "display (full) name of generated library and exit. "
            "MODE = 'load' displays the full file name of the library. "
            "MODE = 'build' displays the name of the library to be "
            "linked to a program using that library. "
        )
    return parser


COMPILER_ERROR = "Don't know how to deal with %s compiler"


def get_environment(name):
    data = []
    try:
        value = os.environ[name]
        values = value.split()
        for s in values:
            s = s.strip()
            if s.isspace():
                continue
            if s[0] == s[-1] == '"':
                s = s[1:-1]  
            data.append(s)  
    except KeyError:
        pass
    return data


def c_define_args(cmdline_args):
    compiler = cmdline_args.compiler
    cargs = []
    if compiler in ['mingw32', 'unix', 'darwin']:
        def_, undef_ =  '-D', '-U'
    elif compiler in ['msvc']:
        def_, undef_ =  '/D', '/U'
    else:
        raise ValueError(COMPILER_ERROR % compiler)  
    # add define macros
    for name in cmdline_args.define:
        cargs.append(def_ + name)
    # add undef macros
    for name in cmdline_args.undef:
        cargs.append(undef_ + name)
    return cargs


def linker_library_args(cmdline_args):
    compiler = cmdline_args.compiler
    if compiler in ['unix', 'mingw32', 'darwin']:
        path, lib = "-L", "-l"
    elif compiler in ['msvc']:
        path, lib = "/LIBPATH:", "/DEFAULTLIB:"
    else:
        raise ValueError(COMPILER_ERROR % compiler)
    largs = []
    for dir in cmdline_args.library_path:
        largs.append(path + dir)
    if compiler in ['mingw32']:
        # special hack between unix and windows world
        for library in cmdline_args.libraries:
            if library.startswith('lib'):
                library = library[3:]
            if library.endswith('.lib'):
                library = library[:-4]
            largs.append(lib + library )
    else:
        for library in cmdline_args.libraries: 
            largs.append(lib + library )
    return largs


def path_from_arg(cmdline_args, path_arg, filename, force_path = False, flat = False):
    path = getattr(cmdline_args, path_arg, None)
    if not path:
        path = ""
        if force_path:
            ERR = "Missing command line argument '--%s'"
            raise ValueError(ERR % re.sub('_', '-', path_arg))
    else:        
        if force_path:
            real_path = os.path.realpath(os.path.normpath(path))
            os.makedirs(real_path, exist_ok=True)
    arg_list = os.path.split(filename)
    if len(arg_list) < 1:
        ERR = "Bad file name for command line argument '--%s'"
        raise ValueError(ERR % re.sub('_', '-', path_arg))
    if flat:
        arg_list = arg_list[-1:]
    new_path = os.path.join(path, *arg_list)
    return os.path.realpath(os.path.normpath(new_path))
     


def make_source_object_pairs(cmdline_args):
    if os.name == "posix":
        ext = '.o'
    elif os.name == "nt":
        ext = '.obj'
    else:
        ERR = "Dont't know object file extension in platform %s"
        raise ValueError(ERR % os.name)
    output = []
    for source in cmdline_args.sources:
        c_name = path_from_arg(cmdline_args, 'source_dir', source)
        objname = os.path.splitext(os.path.split(source)[-1])[0] + ext
        o_name = path_from_arg(cmdline_args, 'library_dir', objname,
              force_path=True, flat=True)
        output.append([c_name, o_name])
    return output


def process_flags(arg):
    args =  [x.strip() for x in arg.split(',')]
    return [x for x in args if len(x)]


def output_names(cmdline_args):
    name, static = cmdline_args.name, cmdline_args.static
    path_arg = 'library_dir' if static else 'shared_dir'
    lib_name = shared_lib_name(name, 'load', static=static, flat=True)
    lib_path = path_from_arg(cmdline_args, path_arg, lib_name,
         force_path = True)
    if os.name == "nt" and not static:
        implib_name = shared_lib_name(name, 'build', static=static, 
             flat=True)
        implib_path = path_from_arg(cmdline_args, 'library_dir', 
             implib_name, force_path = True)
    else:
        implib_path = None
    return lib_path, implib_path
    

def set_ld_library_path_args(cmdline):
    args = []
    for path in cmdline.rpath:
        args.append("-Wl,-rpath," + path)
    return args

def make_dll_nt_msvc(cmdline_args):
    """Create a Windows DLL with the mingw compiler"""
    compile_args = ["cl", "/c", "/O2",  "/W4", "/DMS_WIN64"]
    compile_args += process_flags(cmdline_args.cflags)
    for ipath in cmdline_args.include_path:
        compile_args.append("/I" + os.path.realpath(ipath))
    compile_args += c_define_args(cmdline_args)
    objects = []
    arglist = []
    # Compile sources and add objects to list 'objects'
    for source, obj in make_source_object_pairs(cmdline_args):
        args = compile_args + [source, '/Fo%s' % obj]
        arglist.append(args)
        objects.append(obj)
    workers = SimpleProcessWorkers(cmdline_args.n)
    workers.run(arglist)

    # Link
    lib, implib = output_names(cmdline_args)
    if cmdline_args.static:
        lcmd =  ["lib"] + objects + ["/OUT:" + lib ]
    else:
        lcmd = ["link", "/DLL"] + objects
        lcmd += process_flags(cmdline_args.lflags)
        lcmd += linker_library_args(cmdline_args)
        lcmd += ["/OUT:" + lib ]
        lcmd += [ "/IMPLIB:" + implib]
    print(" ".join(lcmd))
    subprocess.check_call(" ".join(lcmd)) 
    print(lib + "\nhas been generated.\n")





def make_so_posix_gcc(cmdline_args):
    """Create a posix shared library with gcc"""
    compile_args = ["cc"]
    compile_args += get_environment("CFLAGS")
    compile_args += ["-c", "-O3", "-Wall"]
    compile_args += process_flags(cmdline_args.cflags)
    for ipath in cmdline_args.include_path:
        compile_args += ["-I", os.path.realpath(ipath)]
    compile_args += c_define_args(cmdline_args)
    objects = []
    arglist = []
    # Compile sources and add objects to list 'objects'
    for source, obj in make_source_object_pairs(cmdline_args):
        args = compile_args + ["-fPIC", source, "-o", obj]
        arglist.append(args)
        objects.append(obj)
    workers = SimpleProcessWorkers(cmdline_args.n)
    workers.run(arglist)

    # Link
    lib, implib = output_names(cmdline_args)
    if cmdline_args.static:
        lcmd =  ["ar"]
        lcmd +=  ["rcs", lib ] + objects
    else:
        lcmd = ["cc"]
        lcmd += get_environment("LDFLAGS")
        lcmd += ["-shared",  "-Wall"]
        lcmd += process_flags(cmdline_args.lflags)
        lcmd += objects + linker_library_args(cmdline_args)
        lcmd += ["-o", lib ]
    lcmd += set_ld_library_path_args(cmdline_args)
    print(" ".join(lcmd))
    subprocess.check_call(lcmd) 
    print(lib + "\nhas been generated.\n")


def make_dll_nt_mingw32(cmdline_args):
    """Create a Windows DLL with the mingw compiler"""
    compile_args = ["gcc"]
    compile_args += get_environment("CFLAGS")
    compile_args += ["-c", "-O3", "-Wall", "-DMS_WIN64"]
    if cmdline_args.static:
         compile_args.append("-mno-stack-arg-probe")
    # Option  "-mno-stack-arg-probe" prevents the linker error
    # that the symbol ___chkstk_ms  cannot be found.
    compile_args += process_flags(cmdline_args.cflags)
    for ipath in cmdline_args.include_path:
        compile_args += ["-I", os.path.realpath(ipath)]
    compile_args += c_define_args(cmdline_args)
    objects = []
    arglist = []
    # Compile sources and add objects to list 'objects'
    for source, obj in make_source_object_pairs(cmdline_args):
        args = compile_args + [source, '-o', obj]
        arglist.append(args)
        objects.append(obj)
    workers = SimpleProcessWorkers(cmdline_args.n)
    workers.run(arglist)

    # Link
    lib, implib = output_names(cmdline_args)
    if cmdline_args.static:
        lcmd =  ["ar", "rcs", lib ] + objects
    else:
        lcmd = ["gcc"] 
        lcmd += get_environment("LDFLAGS")
        lcmd += objects
        lcmd += process_flags(cmdline_args.lflags)
        lcmd += linker_library_args(cmdline_args)
        lcmd +=  ["-o",  lib, "-s", "-shared"]
        lcmd += ["-Wl,--subsystem,windows,--out-implib," + implib]
    print(" ".join(lcmd))
    subprocess.check_call(lcmd) 
    print(lib + "\nhas been generated.\n")


def make_so_darwin_gcc(cmdline_args):
    """Create a macOS shared library with gcc"""
    return make_so_posix_gcc(cmdline_args)



def build_shared_library(cmdline_args):
    """Build a shared library

    """ 
    ERR_BLD = "Don't know how to build a %s with the '%s' compiler"
    if cmdline_args.mockup:
        return
    compiler = cmdline_args.compiler
    if os.name == "posix":
        if compiler == 'unix':
            make_so_posix_gcc(cmdline_args)
        elif compiler == 'darwin':
            make_so_darwin_gcc(cmdline_args)
        else:
            raise ValueError(ERR_BLD % ("posix shared library", compiler))
    elif os.name == "nt":
        if compiler == 'mingw32':
            make_dll_nt_mingw32(cmdline_args)
        elif compiler == 'msvc':
            make_dll_nt_msvc(cmdline_args)
        else:
            raise ValueError(ERR_BLD % ("Windows DLL", compiler))
    else:
        ERR = "Don't know how to build a shared library on platform '%s'" 
        raise ValueError(ERR % os.name)




def build(args):
    parser = build_shared_lib_parser()
    cmdline_args = parser.parse_args(args)
    build_shared_library(cmdline_args)


if __name__ == "__main__":
    build(sys.argv[1:])


  
