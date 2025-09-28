"""Auxiliary functions for writing tables into C program"""

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import types
import sys
import re
import os
import shutil
import glob
import warnings
from io import IOBase
import argparse
import textwrap
from collections import OrderedDict
from collections.abc import Mapping
from importlib import import_module
import warnings

from mmgroup.generate_c.make_c_tables import TableGenerator
from mmgroup.generate_c.make_c_tables import NoDirectives


class GenerateAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if getattr(namespace, self.dest, None) is None:
            setattr(namespace, self.dest, [])
        d = getattr(namespace, self.dest, None)
        d.append((option_string[2:], values))

class GenerateTables(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if getattr(namespace, self.dest, None) is None:
            setattr(namespace, self.dest, [])
        d = getattr(namespace, self.dest, None)
        param_lists = getattr(namespace, 'param', None)
        params = param_lists[-1] if param_lists else []
        for v in values:
            d.append((v, params))

class MyArgumentParser(argparse.ArgumentParser):
    def convert_arg_line_to_args(self, arg_line):
        return arg_line.split()



def generate_code_parser():
    description = ('Generate C files and a header file for the mmgroup project '
    'from source files using certain tables and directives in the source file. '
    'Optionally, a .pxd, and a .pxi file may be generated from that header. ' 
    )

    # epilog = ("Some more documentation will follow here later."
    #)

    parser = MyArgumentParser(
        prog = 'generate_code',
        description=description,  
        # epilog=epilog,
        fromfile_prefix_chars='+',
    )

    parser.add_argument('--sources', 
        dest='actions', nargs = '*',  metavar='SOURCE',
        action=GenerateAction,
        help="List of SOURCE files to be generated. For each SOURCE with "
             "extension '.c' a '.c' file is generated from a file with "
             "the same name and extension '.ske'. A SOURCE with "
             "extension '.h' is copied into the common header file. "
             "Each SOURCE is seached in the path set by parameter "
             "'--source-path'. Output is written to the directory set "
             "by parameter '--out-dir'."
    )

    parser.add_argument('--out-header', 
        metavar = 'HEADER', default = None, nargs = '*',
        help = "Set name of output header file to HEADER. By default we take "
               "the first file with extension .h in the list given in the "
               "argument '--sources'." )

    parser.add_argument('--copy',
        nargs = '*',  metavar='FILE',
        action = 'extend', default = [],
        help = "Copy FILE to directory set by parameter '--out-dir'. "
               "Each FILE is searched in the path set by parameter "
               "'--source-path'. Wildcards in FILE are allowed."
    )

    parser.add_argument('--pxd', 
        metavar='PXD',
        action='store', default = None,
        help = "Set input '.pxd' file PXD for generating '.pxd' file with "
               "same name from that input file and from the generated header."
    )

    parser.add_argument('--pxi', 
        action='store_true', 
        help="Generate '.pxi' file from '.pxd' file if this option is set."
    )

    parser.add_argument('--pyx', 
        metavar='PYX',
        action='store', default = None,
        help="Copy input '.pyx' file PYX from source path to output directory."
    )

 
    parser.add_argument('--tables', 
        nargs = '*',  metavar='TABLES',
        action = 'extend', default = [], 
        help="Add list TABLES of table classes "
        "to the tables to be used for code generation."
    )


    parser.add_argument('--set', 
        dest='actions',
        nargs = '+',  metavar='VAR=VALUE',
        action=GenerateAction,
        help="Set variable VAR to value VALUE. "
        "When generating code with subsequent '--sources' options "
        "then the table classes will set VAR=VALUE."
    )

    parser.add_argument('--subst', 
        dest='actions',
        nargs = '*',  metavar='PATTERN SUBSTITUTION',
        action=GenerateAction,
        help=textwrap.dedent(
        "Map the name of a '.c' or '.h' file to be generated to the "
        "name of a file used as a source for the generation process. "
        "\n"
        "Example: \"--subst mm(?P<p>[0-9]+)_op mm_op\"   maps e.g."
        "'mm3_op' to 'mm_op'. "
        "We substitute PATTERN in the name of a generated file by "
        "SUBSTITUTION for obtaining the name of that source. "
        "The part \"(?P<p>[0-9]+)\" "
        "creates a variable 'p' that takes the decimal string "
        "following the initial letters 'mm' in the file name. "
        "Then variable 'p' will be passed to the table classes "
        "used for code generation. "
        "PATTERN must be given in python regular expression syntax. "
        )
    )
   
    parser.add_argument('--source-path', nargs = '*', action='extend',
        metavar = 'PATHS', default = [],
        help = 'Set list of PATHS for finding source files to be used '
               'for generation process.')

    parser.add_argument('--py-path', nargs = '*', action='extend',
        metavar = 'PATHS', default = [],
        help = 'Set list of PATHS for finding python scripts')

    parser.add_argument('--out-dir', 
        metavar = 'DIR', default = None,  action = 'store',
        help = 'Store output files with extensions .c and .h '
       'in directory DIR.')

    parser.add_argument('--out-pxd-dir', 
        metavar = 'DIR', action='store', default = None,
        help = 'Store output files with extensions .pxd, .pxi, and .pyx '
       'in directory DIR.')

    parser.add_argument('--dll', 
        default = None,
        action = 'store',
        metavar = 'DLL_NAME',
        help = 'Generate code for exporting C functions to a DLL '
               'or to a shared library with name DLL_NAME. Parameter '
               'DLL_NAME must be the same for all generated C files '
               'to be placed into the same DLL. '
               'This parameter is not used for any other purposes.'
               '\'--dll None\' generates code for static linking.' )

    parser.add_argument('--nogil', 
        action='store_true', default = False,
        help = "Optional, declare functions from .pxi file as 'nogil' "
        "when set."
    )

    #parser.add_argument('--n',
    #    type=int,  default = 1, metavar = 'N',
    #    help = 'Use N parallel processes for compiling.'
    #    )

    parser.add_argument('--mockup', 
        default = False,
        action = 'store_true',
        help = 'Use tables and directives for Sphinx mockup if present')
 

    parser.add_argument('--library-path', nargs = '*', action='extend',
        metavar = 'PATHS', default = [],
        help = 'Set list of PATHS for finding shared libraries')

    parser.add_argument('--no-library-path', action = 'store_true',
        help=argparse.SUPPRESS)

    parser.add_argument('-v', '--verbose', action = 'store_true',
        help = 'Verbose operation')

    return parser



def set_real_path(instance, attribute):
    path = getattr(instance, attribute)
    if path is None:
        return
    real_path = os.path.realpath(path)
    setattr(instance, attribute, real_path) 

def set_real_pathlist(instance, attribute):
    pathlist = getattr(instance, attribute)
    if pathlist is None:
        return
    if not isinstance(pathlist, list):
         ERR = "Attribute '%s' of class %s object is not a list of paths"
         raise TypeError(ERR % (attribute, instance.__class__))
    real_pathlist = [os.path.realpath(path) for path in pathlist]
    setattr(instance, attribute, real_pathlist) 

def find_file(pathlist, filename, verbose = 0):
    if verbose:
        print("trying to find file", filename)
        print("pathlist is", pathlist)
    fpath, fname = os.path.split(filename)
    fpath = os.path.normpath(fpath) if len(fpath) else ""
    filename = os.path.join(fpath, fname)
    for path in pathlist:
        try:
            if verbose:
                print("inspecting", os.path.join(path, fpath))
            files = os.listdir(os.path.join(path, fpath))
            if verbose:
                print("files are", files)
            if fname in files:
                if verbose:
                    print("found", os.path.join(path, filename)) 
                return os.path.join(path, filename)
        except:
            pass
    ERR = "File %s not found" 
    raise IOError(ERR % filename)
     
def set_real_out_header(instance):   
    dir = getattr(instance,'out_dir')
    header =  getattr(instance,'out_header')
    if header is None:
        setattr(instance, 'real_out_header', None)
        return
    real_header = os.path.join(dir, os.path.normpath(header))
    setattr(instance, 'real_out_header', real_header) 
     

def search_files(patterns, directories):
    """Search files matching given wildcard patterns in given directories.

    Args:
        patterns (list): List of wildcard patterns (e.g., ['*.txt', '*.py'])
        directories (list): List of directories to search in.

    Returns:
        list: List of tuples (directory, relative file path).
    """
    result = []
    for pattern in patterns:
        for directory in directories:
            dir = os.path.normpath(directory)
            search_path = os.path.join(dir, pattern)
            for match in glob.glob(search_path, recursive=True):
                relative_path = os.path.relpath(match, dir)
                result.append((dir, relative_path))
    return result


def copy_files_with_dirs(src_pattern, src_dirs, dest_dir, verbose = 0):
    """ Copy a source files to a destination path.

    Args:
        src_pattern (list): List of wildcard patterns for source files.
        src_dir (list): List of directories to search for source files.
        dest_dir (str): directory where to store the copied files
    """
    result = search_files(src_pattern, src_dirs)
    dest_dir =  os.path.normpath(dest_dir)
    for src_dir, relative_path in  result:
        src = os.path.join(src_dir, relative_path)
        dest = os.path.join(dest_dir, relative_path)
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        shutil.copy2(src, dest)
        if verbose:
           print(f"Copied '{src}' to '{dest}'")





class ActivatedPythonPath:
    """Add a list of paths to the path where python finds scripts

    The constructor of this class takes a list of names of paths
    to be added to the python path given by ``sys.path``. Calling

    ``pth = ActivatedPythonPath(paths)``

    adds the paths in the list ``path`` to ``sys.path``.

    You may call ``pth.close()`` for reversing this operation.
    """
    old_py_path = []
    new_py_path = []

    def __init__(self, paths = None):
        def purge_paths(paths):
            if not paths:
                return []
            if isinstance(paths, str):
                return [paths]
            for s in paths:
                assert isinstance(s, str)
            return paths
        self.activate_py_path(purge_paths(paths))


    def activate_py_path(self, directories):
        self.old_py_path = [x for x in sys.path]
        self.new_py_path = []
        for i, path in enumerate(directories):
            sys.path.insert(i, path)
            self.new_py_path.append(path)

    def deactivate_py_path(self):
        if len(self.new_py_path) == 0:
            return
        if sys.path != self.new_py_path + self.old_path:
            W = ("Could not deactivate python path activated " 
                 "by class ActivatedPythonPath object")
            warnings.warn(W, UserWarning)
        else:
            del sys.path[:len(self.s.py_path)]

    def close(self):
        self.de_activate_py_path()


class ImmutableDict(Mapping):
    """Create an immutable and hashable snapshot of a dictionary

    """
    def __init__(self, *dictionaries):
        d = {}
        for dictionary in dictionaries:
            if dictionary:
                d.update(dictionary)
        self.d = OrderedDict()
        for key in sorted(d.keys()):
            self.d[key] = d[key] 
        self._hash =  hash(tuple(self.d.items()))
    def __hash__(self):
        return self._hash  
    def __getitem__(self, key):
        return self.d[key]
    def __iter__(self):
        yield from self.d
    def __len__(self):
        return len(self.d)
    def str(self):
        return "Immutable" + str(dict(self.d))
    __repr__ = str

EmptyImmutableDict = ImmutableDict({})


ERR_SET = "Arguments following option --set must have shape VAR=VALUE" 
ERR_SUBST = "Two arguments PATTERN, SUBSTITUTION must follow option --subst" 

m_set = re.compile(r"([A-Za-z][A-Za-z0-9_]*)=(.*)")



"""
def subst_C_name(name, subst, param, extension):
    dir, filename = os.path.split(name)
    name, _extension = os.path.splitext(filename)
    if _extension == ".h":
        return None
    if subst:
        #print("substitute", name, "using", subst, param) 
        name = re.sub(subst[0], subst[1], name)
        if param:
            name = name.format(**param)
    return os.path.join(dir, name + '.' + extension)    
"""

def subst_source_name(input_name, subst):
    dir, filename = os.path.split(os.path.normpath(input_name))
    name, ext = os.path.splitext(filename)
    #print(dir, name, ext)
    if ext == '.c':
        source_ext = '.ske'
        target_name = name
    elif ext == '.h':
        source_ext = '.h'
        target_name = None
    else:
        ERR = "Don't know how to generate source from file"
        raise ValueError(ERR + filename)
    mm = re.compile(subst[0]) if subst and len(subst) else None
    if mm and len(subst) > 1:
        source_name = mm.sub(subst[1], name)
    else:
        source_name = name
    m = mm.match(name) if mm else None
    vars = m.groupdict() if m else {}
    source = os.path.join(dir, source_name + source_ext)
    source = os.path.normpath(source)
    if target_name:
        target = os.path.join(dir, target_name + ext)
        target = os.path.normpath(target)
    else:
        target = None
    return source, target, ext, vars    


def make_actions(s):
    n = 1
    param = {}
    s.c_files = OrderedDict()
    subst = None
    if not s.actions:
        return
    for action, data in s.actions:
        if action == 'set':
            for instruction in data:
                m = m_set.match(instruction)
                if m is None:
                    raise ValueError(ERR_SET)
                var, value = m.groups()
                value = value.strip()
                if len(value):
                    param[var] = value
                else:
                    del param[var]
        if action == 'sources':
            for name in data:
                src, target, ext, vars = subst_source_name(name, subst)
                if not s.out_header and ext == '.h':
                    s.out_header = name 
                param_dict = ImmutableDict(param, vars)
                if param_dict not in s.c_files:
                    s.c_files[param_dict] = []
                s.c_files[param_dict].append((src, target, n))
                n += 1
        if action == 'subst':
            if len(data) == 0:
                subst = None   
            elif len(data) != 2:
                raise ValueError(ERR_SUBST)
            subst = data
 
            
def finalize_parse_args(s):
    #print("\nFinalizing\n", s)
    if getattr(s, 'finalized', None):
        return
    set_real_pathlist(s, 'source_path')
    set_real_pathlist(s, 'py_path')
    set_real_path(s, 'out_dir')
    make_actions(s)
    set_real_out_header(s)
    #print("\nFinalized\n", s)
    s.finalized = True
 

def import_tables(table_modules, mockup = False, verbose = False):
    """Import tables from python modules
    
    Here ``table_modules`` is a list of names of python modules in
    the usual python syntax. This function tries to import a class 
    with name ``Tables`` from each of these modules. 
    It returns the list of imported classes.

    If ``mockup`` is True then the function tries to import a class 
    with name ``MockupTables`` instead of class ``Tables``. If this
    fails then it tries to import class  ``Tables``. 

    """
    table_classes =  []
    if table_modules:
        for module in table_modules:
            if verbose:
                print("Importing Tables from module", module)
            m = import_module(module)
            if mockup:
                try:
                    table_class = m.MockupTables
                except:
                    table_class = m.Tables
            else:
                table_class = m.Tables
            table_classes.append(table_class)
    return table_classes    


def load_tables(tg, tables, params, directives=True):
    """Load tables into instance ``tg`` of class ``TableGenerator``

    The  argument  ``tables`` must be a list of table_classes.
    Such a table class should be a class with attributes ``tables``
    and ``directives`` that are dictionaries mapping names to tables
    or to directives. 

    Argument ``params`` should be a mapping from parameter
    names to values.

    Then for each table class the parameters occuring in ``params``
    are set to the values given by the mapping ``params``; and the 
    dictionary mapping these parameters to their values is passed 
    to the constructor of class ``table_class`` as keyword arguments.
    We recommend that the constructor of a table class accepts
    arbitrary keyword arguments and ignores redundant arguments.

    The instance of ``table_class`` constructed that way should also 
    have dictionaries ``tables`` and ``directives`` as attributes as
    above. The tables and directives of instance ``tg`` of class
    ``TableGenerator`` are initialized and the updated with all the
    corresponding dictionaries. Afterwards ``tg`` may be used for code
    generation

    If the arguments ``directives`` is False then the table
    generator ``tg`` will support no directives at all.  
    """
    assert isinstance(params, Mapping)
    assert isinstance(tg, TableGenerator)
    _tables = {}
    _directives = {} if directives else NoDirectives
    for table_class in tables:
        m_tables = table_class(**params)
        #print("Loading", table_class, params)
        try:
            new_tables = m_tables.tables
            _tables.update(new_tables)
        except:
            print("Could not load table from class", table_class)
            print("Parameters:", dict(params))
            new_tables = getattr(m_tables, 'tables', "<undefined>")
            print("Attribute 'tables' of table class is:", new_tables)
            raise
        if directives:
            new_directives = m_tables.directives
            _directives.update(new_directives)
    tg.set_tables(_tables, _directives)


m_split_kwd =  re.compile(r"\s*//\s*\%\%INCLUDE_HEADERS")


class StringOutputFile():
    """Simulate a text stream open for output.

    Suports method ``write`` for text streams only, writing
    data to a internal buffer.

    After writing, the stream can be read by standard iterator
    methods in the same way as a file being opened for reading.
    """
    def __init__(self):
        self.s = []
        self.output_closed = False
        self.pos = 0
    def write(self, data):
        if self.output_closed:
           ERR = "Attempt to write to StringOutputFile after reading"
           raise IOError(ERR)
        if data:
            for line in data.splitlines(keepends=True):
                 self.s.append(line)
    def copy_to(self, ostream):
        """Copy internal data to test_stream ``iostream``"""
        self.output_closed = True
        for s in self.s:
             ostream.write(s)
    def as_string(self):
        self.output_closed = True
        return "".join(self.s)    
    def close(self):
        self.output_closed = True
    def __iter__(self):
        self.output_closed = True
        return self
    def __next__(self):
        assert self.output_closed
        if self.pos >= len(self.s):
           raise StopIteration
        self.pos += 1
        return self.s[self.pos - 1]

class NullOutputFile:
    """Simulate a forgetful text stream open for output."""
    def write(self, *args, **kwds):
        return
    close = write


def split_header(file):
    head = StringOutputFile()
    tail = StringOutputFile()
    split_comment = ""
    if not file:
         return head, split_comment, tail
    for line in file:
        if m_split_kwd.match(line):
             split_comment = line
             break
        else:
             head.write(line)
    for line in file:
        tail.write(line)
    head.close()
    tail.close()
    return head, split_comment, tail


def open_for_write(dir, filename):
    """Open file ``dir\filename`` for output""" 
    if not filename:
        return NullOutputFile
    pathname = os.path.join(dir, filename)
    dir_name = os.path.split(pathname)[0]
    if len(dir_name):
        os.makedirs(dir_name, exist_ok=True)
    return open(pathname, "wt")

def write_warning_generated(stream):
    stream.write(
    "// Warning: This file has been generated automatically."
    " Do not change!\n" 
    )


def check_arg_parser(parser_namespace):
    p = parser_namespace
    ok = False
    try:
        ok = (type(p.actions) == list and
            type(p.tables) == list)
    except:
        pass
    if not ok:
        ERR = "Error in argument parser for code generation" 
        raise ValueError(ERR)


############################################################################
# Dealing with DLLs
############################################################################

DLL_HEADER_PREFIX = """

/// @cond DO_NOT_DOCUMENT 
//  Definitions for using this header in a a DLL (or a shared library)

// Generic helper definitions for DLL (or shared library) support
#if defined(_WIN32) || defined(__CYGWIN__)
  #define {0}_DLL_IMPORT __declspec(dllimport)
  #define {0}_DLL_EXPORT __declspec(dllexport)
#elif (defined(__GNUC__) || defined(__clang__)) && defined(_WIN32)
  #define {0}_DLL_IMPORT __attribute__((noinline,optimize("no-tree-vectorize"),visiblity("default")))
  #define {0}_DLL_EXPORT __attribute__((noinline,optimize("no-tree-vectorize"),visiblity("default")))
#else
  #define {0}_DLL_IMPORT
  #define {0}_DLL_EXPORT
#endif

// Now we use the generic helper definitions above to define {0}_API 
// {0}_API is used for the public API symbols. It either DLL imports 
// or DLL exports 

#ifdef {0}_DLL_EXPORTS // defined if we are building the {0} DLL 
  #define {0}_API {0}_DLL_EXPORT
#else                  // not defined if we are using the {0} DLL 
  #define {0}_API  {0}_DLL_IMPORT
#endif // {0}_DLL_EXPORTS

/// @endcond
"""



DLL_C_FILE_PREFIX = """
/// @cond DO_NOT_DOCUMENT 
#define {0}_DLL_EXPORTS 
/// @endcond

"""

EXPORT_KWD = "{0}_API" 

############################################################################
# class CodeGenerator
############################################################################


class CodeGenerator:
    """Yet to be documented"""
    activated_python_path = None
    BIGINT = 0x7fffffff
    def __init__(self, args):
        if isinstance(args, str):
            parser = generate_code_parser()
            self.s = parser.parse_args(args.split())
        elif isinstance(args, list):
            parser = generate_code_parser()
            self.s = parser.parse_args(args)
        elif isinstance(args, argparse.Namespace):
            check_arg_parser(args)
            self.s = args
        else:
            ERR = "Cannot construct class %s object from %s object"
            raise TypeError(ERR % (self.__class__, type(args)))
        self.old_path = None
        finalize_parse_args(self.s)
        self.tg = TableGenerator()

    def check_input_files(self):
        finalize_parse_args(self.s)

    def find_file(self, filename):
        return find_file(self.s.source_path, filename)

    def activate_py_path(self):
        finalize_parse_args(self.s)
        self.activated_python_path = ActivatedPythonPath(self.s.py_path)

    def deactivate_py_path(self):
        if self.activated_python_path is not None:
            self.activated_python_path.close()
            self.activated_python_path = None

    def import_tables(self):
        finalize_parse_args(self.s)
        if getattr(self.s, "table_classes", None) is None:
            table_modules = self.s.tables
            verbose = self.s.verbose
            mockup = self.s.mockup
            table_classes = import_tables(table_modules, mockup, verbose)
            self.s.table_classes = table_classes    
        
    def load_table_generator(self, table_generator, params):
        assert isinstance(params, ImmutableDict)
        assert isinstance(table_generator, TableGenerator)
        #tables = {}
        #directives = {}
        tables = self.s.table_classes
        load_tables(table_generator, tables, params)

    def set_dll_prefixes(self, table_generator):
        if self.s.dll and self.s.dll != "None" and not self.s.mockup:
            self.c_prefix = DLL_C_FILE_PREFIX.format(self.s.dll)
            self.h_prefix = DLL_HEADER_PREFIX.format(self.s.dll)
            table_generator.export_kwd = EXPORT_KWD.format(self.s.dll)
        else:
            self.c_prefix = ""
            self.h_prefix = ""

    def copy_files(self):
        copy_files_with_dirs(self.s.copy, self.s.source_path,
            self.s.out_dir, self.s.verbose)

    def generate_c(self, src_name, dest_name, n, h_prefix=""):
        src_file_name = self.find_file(src_name)
        src = open(src_file_name, "rt")
        s = self.s
        if dest_name:
            dest = open_for_write(s.out_dir, dest_name)
            if s.verbose:
                 print("Generating file %s" % dest_name)
            write_warning_generated(dest)
            out_h = StringOutputFile()
            out_h.write(h_prefix)
            out_h.write("// %%%%FROM %s\n" % dest_name)
            dest.write(self.c_prefix)
            self.tg.generate(src, dest, out_h)
            out_h.write("\n")
            dest.close()
            log = "Inserting header for " + dest_name
            return [(n, out_h.as_string(), log)]
        else:
            out_h = StringOutputFile()
            end_h = StringOutputFile()
            h_head, h_split, h_tail = split_header(src)
            self.tg.generate(h_head, None, out_h, 'h')
            out_h.write("\n")
            end_h.write(h_split)
            self.tg.generate(h_tail, None, end_h, 'h')
            end_h.write("\n")
            log = "Inserting header " + src_file_name
            return [
               (n, out_h.as_string(), log),
               (self.BIGINT-n, end_h.as_string(), "")
            ]

    def generate_c_files(self, gen_list):
        h_list = []
        for args in gen_list:
            h_list += self.generate_c(*args)
        return h_list
        """ # The following does work due to pickle problems
        h_list = []
        processes = self.s.n
        if processes > 1:
            from multiprocessing import Pool, cpu_count
            if cpu_count():
                processes = min(processes, cpu_count())
        if processes <= 1:
            for args in gen_list:
                h_list += self.generate_c(*args)
            return h_list
        with Pool(processes = processes) as pool:
             results = pool.starmap(self.generate_c, gen_list)
        pool.join()
        for data in results:
            h_list += data
        return h_list
        """

    def generate(self):
        finalize_parse_args(self.s)
        h_list = []
        s = self.s
        out_dir = s.out_dir
        self.set_dll_prefixes(self.tg)
        if s.verbose:
            print("Generating header %s" % s.out_header)
        for param, c_files in s.c_files.items():
            gen_list = []
            self.load_table_generator(self.tg, param)
            for src_name, dest_name, n in c_files:
                h_prefix = self.h_prefix if dest_name else ""
                gen_list.append((src_name, dest_name, n, h_prefix))
                self.h_prefix = "" if  h_prefix else self.h_prefix
            h_list += self.generate_c_files(gen_list)
        if s.real_out_header:
            out_header = open_for_write(out_dir, s.real_out_header)
            write_warning_generated(out_header)
            for i, text, log in sorted(h_list):
                if s.verbose and log:
                      print(log)
                out_header.write(text)
            out_header.close()

        self.copy_files()
 
        if not s.mockup:
            from mmgroup.generate_c.generate_pxd import make_pxd
            self.load_table_generator(self.tg, ImmutableDict())
            make_pxd(s, self.tg)

    def c_files(self):
        finalize_parse_args(self.s)
        result = []
        for _, c_files in self.s.c_files.items():
            for _, dest_name, _1 in c_files:
                if dest_name:
                    unix_dest_name = re.sub(r"\\", "/", dest_name)
                    result.append(unix_dest_name)
        return result
       
            
    def display_args(self):
        s = self.s
        actions = getattr(s, 'actions', None)
        if actions:
            print("actions:")
            for name, action in actions:
                print(" ", name, action)
        c_files = getattr(s, 'c_files', None)
        if c_files:
            for par, files in c_files.items():
                dict = {}
                dict.update(par)
                print("c_files (%s):" % dict)
                for src, dest, n in files:
                    print(" ", src, ",", dest, "," , n)
        tables = getattr(s, 'tables', None)
        if tables:
            print("tables:")
            for table in tables:
                print(" ", table)

        attribs = [
            'source_header', 'out_header', 'py_path', 'library_path',
            'out_dir', 'verbose', 'dll'
        ]
        for attr in attribs:
            print(attr + ':', getattr(s, attr, '<undefined>'))

############################################################################
# Auxiliary functions
############################################################################
           

def set_shared_libraries(parsed_args):
    """Set path for finding shared libraries

    Parameter 'parsed_args' should be an object returned by method

         generate_code_parser().parse_args(args)

    where 'args' is an string containing an argument list.

    The function extends environment variable 'LD_IBRARY_PATH' (in a
    posix system) or 'path' (in a Windows system) according to that
    argument list, so that shared or dynamic libraries given in that
    list can be found.

    The function returns ``True`` if it has changed the environment.
    If this is the case then the calling function should launch a
    subprocess for further actions, since changing the enviroment
    affects subprocesses of the calling process only.
    """
    ld_args = getattr(parsed_args, 'library_path', [])
    if getattr(parsed_args, 'no_library_path', 0) or len(ld_args) == 0:
        return False;
    if os.name == 'posix':
        ld_path = os.getenv('LD_LIBRARY_PATH')
        ld_path = [] if ld_path is None else [ld_path]
        os.environ['LD_LIBRARY_PATH'] = ':'.join(ld_path + ld_args)
        return True
    elif os.name == 'nt':
        ld_path = os.getenv('path')
        ld_path = [] if ld_path is None else [ld_path]
        os.environ['path'] = ';'.join(ld_path + ld_args)
        return True
    else:
        ERR = "Don't now how to set library path in %s operating system"
        raise ValueError(ERR % os.name)





############################################################################
# Setting shared library path
############################################################################




def parse_set_shared_libraries(args):
    """Set path for finding shared libraries

    Here parameter 'args' should be ``sys.arv[1:]``

    Then the function parses 'args' and performs the same action and
    returns the same value as function ``set_shared_libraries(args)``.
    """
    path_parser = argparse.ArgumentParser(add_help=False)
    path_parser.add_argument('--library-path',
        nargs = '*', action='extend',  default = []
    )
    path_parser.add_argument('--no-library-path', action = 'store_true')
    path_args = path_parser.parse_known_args(args)
    return set_shared_libraries(path_args[0])

