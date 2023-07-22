"""Auxiliary functions for writing tables into C program"""

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import types
import sys
import re
import os
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
    description = ("Generate C code and header files for the mmgroup project "
    "from source files using certain tables and directives in the source file."
    )

    epilog = ("Some more documentation will follow here later."
    )

    parser = MyArgumentParser(
        description=description, epilog=epilog,
        fromfile_prefix_chars='+',
    )

    parser.add_argument('--sources', 
        dest='actions', nargs = '*',  metavar='SOURCES',
        action=GenerateAction,
        help="Add list of SOURCES to actions for generating code. "
             "If one of these SOURCES in the list has extension '.h' "
             "then no corresponding C fle is generated."
    )

    parser.add_argument('--out-header', 
        metavar = 'HEADER', default = None,
        help = 'Set name of output header file to HEADER.')
 
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
        "Example: \"--set P 3  --subst op op{P}\"   maps e.g."
        "'mm_op.c' to 'mm_op3.c'. "
        "Here we substitute PATTERN in the name of a source file by "
        "SUBSTITUTION for obtaining a template for the name of the "
        "corresponding C file. Then we format that template with "
        "Python's .format method (using the values of the current "
        "variables) for obtaining the name of the C file."
        )
    )
   
    parser.add_argument('--source-path', nargs = '*', action='extend',
        metavar = 'PATHS', default = [],
        help = 'Set list of PATHS for finding source files')

    parser.add_argument('--py-path', nargs = '*', action='extend',
        metavar = 'PATHS', default = [],
        help = 'Set list of PATHS for finding python scripts')

    parser.add_argument('--out-dir', 
        metavar = 'DIR', default = None,
        help = 'Set directory DIR for output files')

    parser.add_argument('--export-kwd', 
        default = None,
        action = 'store',
        metavar = 'EXPORT_KWD',
        help = 'Set export keyword for code generator to EXPORT_KWD, deprecated!')


    parser.add_argument('--dll', 
        default = None,
        action = 'store',
        metavar = 'DLL_NAME',
        help = 'Generate code for exporting C functions to a DLL '
               'or to a shared library with name DLL_NAME')

    parser.add_argument('--mockup', 
        default = False,
        action = 'store_true',
        help = 'Use tables and directives for Sphinx mockup if present')
 
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

     
def set_out_header(instance):   
    dir = getattr(instance,'out_dir')
    header =  getattr(instance,'out_header')
    if header is None:
        return
    real_header = os.path.join(dir, os.path.normpath(header))
    setattr(instance, 'out_header', real_header) 
     



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
    def __init__(self, dictionary):
        keys = sorted(dictionary.keys())
        self.d = OrderedDict()
        for key in keys:
            self.d[key] = dictionary[key] 
        self._hash =  hash(tuple(self.d.items()))
    def __hash__(self):
        return self._hash  
    def __getitem__(self, key):
        return self.d[key]
    def __iter__(self):
        yield from self.d
    def __len__(self):
        return len(self.d)
    def restrict(self, keys):
        """Return retriction of this dictionary to a given set of keys""" 
        d = {}
        for key in keys:
            if key in self.d:
                d[key] = self.d[key]
        return d
    def str(self):
        return "Immutable" + str(dict(self.d))
    __repr__ = str

EmptyImmutableDict = ImmutableDict({})


ERR_SET = "Arguments following option --set must have shape VAR=VALUE" 
ERR_SUBST = "Two arguments PATTERN, SUBSTITUTION must follow option --subst" 

m_set = re.compile(r"([A-Za-z][A-Za-z0-9_]*)=(.*)")

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

def make_actions(s):
    n = 1
    param = {}
    s.c_files = OrderedDict()
    subst = None
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
            param_dict = ImmutableDict(param)
            if param_dict not in s.c_files:
                 s.c_files[param_dict] = []
            file_list = s.c_files[param_dict]
            for name in data:
                src_name = os.path.normpath(name)
                dest_name = subst_C_name(name, subst, param, 'c') 
                file_list.append((src_name, dest_name, n))
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
    set_out_header(s)
    make_actions(s)
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
        new_tables = m_tables.tables
        _tables.update(new_tables)
        if directives:
            new_directives = m_tables.directives
            _directives.update(new_directives)
    tg.set_tables(_tables, _directives)


m_split_kwd =  re.compile(r"\s*//\s*\%\%INCLUDE_HEADERS")

def split_header(file):
    head = []
    tail = []
    split_comment = ""
    if not file:
         return head, split_comment, tail
    for line in file:
        if m_split_kwd.match(line):
             split_comment = line
             break
        else:
             head.append(line)
    for line in file:
        tail.append(line)
    return head, split_comment, tail  


class StringOutputFile():
    """Simulate a text stream open for output.

    Suports method ``write`` for text streams only, writing
    data to a internal buffer. 
    """
    def __init__(self):
        self.s = []
    def write(self, data):
        if data:
            self.s.append(data)
    def copy_to(self, ostream):
        """Copy internal data to test_stream ``iostream``"""
        for s in self.s:
             ostream.write(s)
    def as_string(self):
        return "".join(self.s)    
    def close(self):
        self.s = []

class NullOutputFile:
    """Simulate a forgetful text stream open for output."""
    def write(self, *args, **kwds):
        return
    close = write

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
#if defined _WIN32 || defined __CYGWIN__
  #define {0}_DLL_IMPORT __declspec(dllimport)
  #define {0}_DLL_EXPORT __declspec(dllexport)
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
        if self.s.dll and not self.s.mockup:
            self.c_prefix = DLL_C_FILE_PREFIX.format(self.s.dll)
            self.h_prefix = DLL_HEADER_PREFIX.format(self.s.dll)
            table_generator.export_kwd = EXPORT_KWD.format(self.s.dll)
        else:
            self.c_prefix = ""
            self.h_prefix = ""

    def generate(self):
        BIGINT = 0x7fffffff
        finalize_parse_args(self.s)
        out_headers = {}
        s = self.s
        out_dir = s.out_dir
        tg = TableGenerator()
        self.set_dll_prefixes(tg)
        if s.export_kwd:
            #tg.export_kwd = s.export_kwd
            pass
        if s.verbose:
            print("Generating header %s" % s.out_header)
        for param, c_files in s.c_files.items():
            self.load_table_generator(tg, param)
            for src_name, dest_name, n in c_files:
                src_file_name = self.find_file(src_name)
                src = open(src_file_name, "rt")
                if dest_name:
                    dest = open_for_write(out_dir, dest_name)
                    if s.verbose:
                        print("Generating file %s" % dest_name)
                    write_warning_generated(dest)
                    out_h = out_headers[n] = StringOutputFile()
                    out_h.write(self.h_prefix)
                    self.h_prefix = ""
                    out_h.write("// %%%%FROM %s\n" % dest_name)
                    dest.write(self.c_prefix)
                    tg.generate(src, dest, out_h)
                    out_h.write("\n")
                    dest.close()
                else:
                     out_h = out_headers[n] = StringOutputFile()
                     end_h = out_headers[BIGINT-n] = StringOutputFile()
                     h_head, h_split, h_tail = split_header(src)
                     tg.generate(h_head, None, out_h, 'h')
                     out_h.write("\n")
                     end_h.write(h_split)
                     tg.generate(h_tail, None, end_h, 'h')
                     end_h.write("\n")


        out_header = open_for_write(out_dir, s.out_header)
        write_warning_generated(out_header)
        for key in sorted(out_headers):
            out_headers[key].copy_to(out_header)
        out_header.close()

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
            'out_dir', 'export_kwd', 'verbose', 'dll'
        ]
        for attr in attribs:
            print(attr + ':', getattr(s, attr, '<undefined>'))

           


