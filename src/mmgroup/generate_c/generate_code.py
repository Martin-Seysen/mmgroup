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
        help="Add list of SOURCES to actions for generating code."
    )

    parser.add_argument('--source-header', 
       metavar = 'HEADER', # , default = None,
       dest='actions', action=GenerateAction,
       help = 'Set name of source header file to HEADER.')

    parser.add_argument('--out-header', 
        metavar = 'HEADER', default = None,
        help = 'Set name of output header file to HEADER.')
 
    parser.add_argument('--tables', 
        nargs = '*',  metavar='TABLES',
        action=GenerateTables,
        help="Add list TABLES of table classes "
        "to be used for code generation."
    )

    parser.add_argument('--param', 
        nargs = '*',  metavar='PARAM',
        action='append',
        help="Set list PARAM of parameters to be passed to "
        "subsequent table classes. "
        "If no arguments follow then no parameter will be set."
    )

    parser.add_argument('--set', 
        dest='actions',
        nargs = '*',  metavar='PARAM VALUE',
        action=GenerateAction,
        help="Set set parameter PARAM to value VALUE. "
        "When generating code with subsequent '--sources' options "
        "then the table classes will set PARAM=VALUE."
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
        "parameters) for obtaining the name of the C file."
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
        help = 'Set export keyword for code generator to EXPORT_KWD')

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
    raise ValueError(ERR % filename)  

     
def set_out_header(instance):   
    dir = getattr(instance,'out_dir')
    header =  getattr(instance,'out_header')
    if header is None:
        return
    real_header = os.path.join(dir, os.path.normpath(header))
    setattr(instance, 'out_header', real_header) 
     

 



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


ERR_SET = "At most two arguments PARAM, VALUE may follow option --set" 
ERR_SUBST = "Two arguments PATTERN, SUBSTITUTION must follow option --subst" 



def subst_C_name(name, subst, param, extension):
    dir, filename = os.path.split(name)
    name, _extension = os.path.splitext(filename)
    if subst:
        #print("substitute", name, "using", subst, param) 
        name = re.sub(subst[0], subst[1], name)
        if param:
            name = name.format(**param)
    return os.path.join(dir, name + '.' + extension)    

def make_actions(s):
    param = {}
    s.c_files = OrderedDict()
    subst = None
    for action, data in s.actions:
        if action == 'set':
            if len(data) > 2:
                raise ValueError(ERR_SET)
            if len(data) == 2:
                param[data[0]] = data[1]
            if len(data) == 1 and data[0] in param:
                del param[data[0]]
        if action == 'sources':
            file_list = s.c_files.setdefault(ImmutableDict(param), [])
            for name in data:
                src_name = os.path.normpath(name)
                dest_name = subst_C_name(name, subst, param, 'c') 
                file_list.append((src_name, dest_name))
        if action == 'subst':
            if len(data) == 0:
                subst = None   
            elif len(data) != 2:
                raise ValueError(ERR_SUBST)
            subst = data
        if action == 'source-header':
            if len(s.c_files):
                ERR = ( "Option 'source-header' may not occur"
                      "after option 'sources'")
                raise valueError(ERR)
            file_list = s.c_files.setdefault(ImmutableDict(param), [])
            src_name = os.path.normpath(data)
            file_list.append((src_name, None))
 
            


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
 


 

def load_tables(tg, tables, params, mockup=False, directives=True):
    """Load tables into instance ``tg`` of class ``TableGenerator``

    The  argument  ``tables`` must be a list of pairs 
    ``(table_class, table_params)``. The class ``table_class``
    should be a class with attributes ``tables`` and ``directives``
    that are dictionaries mapping names to tables or directives.
    The set ``table_params`` should be a set or a list of strings
    interpreted as names of parameters.

    Argument ``params`` should be a mapping from parameter
    names to values.

    Then for each pair ``(table_class, table_params)`` the parameter
    names occuring in both objects, ``table_params`` and ``params``
    are set to their values given by the mapping ``params``; and 
    a dictionary mapping these parameters to their values is passed
    to the constructor of class ``table_class`` as keyword arguments.

    The instance of ``table_class`` constructed that way should also 
    have dictionaries ``tables`` and ``directives`` as attributes as
    above. The tables and directives of instance ``tg`` of class
    ``TableGenerator`` are initialized and the updated with all the
    corresponding dictionaries. Afterwards ``tg`` may be used for code
    generation

    If ``mockup`` is True the attributes ``tables`` and ``directives``
    of the  instances of  ``table_class`` are replaced by the
    attributes ``mockup_tables`` and ``mockup_directives``, if
    present.

    If the arguments ``directives`` is False then the table
    generator ``tg`` will support no directives at all.  
    """
    assert isinstance(params, Mapping)
    assert isinstance(tg, TableGenerator)
    _tables = {}
    _directives = {} if directives else NoDirectives
    for table_class, table_params in tables:
         common_params = {}
         for param in table_params:
             if param in params:
                 common_params[param] = params[param]
         m_tables = table_class(**common_params)
         if mockup and getattr(m_tables, 'mockup_tables'):
             new_tables = m_tables.mockup_tables
         else:
             new_tables = m_tables.tables
         _tables.update(new_tables)
         if directives:
             if mockup and getattr(m_tables, 'mockup_directives'):
                 new_directives = m_tables.mockup_directives
             else:
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
    return open(os.path.join(dir, filename), "wt")

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

class CodeGenerator:
    """Yet to be documented"""
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
        self.old_path = None
        finalize_parse_args(self.s)

    def check_input_files(self):
        finalize_parse_args(self.s)

    def find_file(self, filename):
        return find_file(self.s.source_path, filename)

    def activate_py_path(self):
        finalize_parse_args(self.s)
        s = self.s
        if not s.py_path:
            return
        if self.old_path:
            ERR = "Python path for class CodeGenerator object already active" 
            raise ValueError(ERR)
        self.old_path = [x for x in sys.path]
        for i, path in enumerate(s.py_path):
            sys.path.insert(i, path)

    def deactivate_py_path(self):
        finalize_parse_args(self.s)
        if not self.s.py_path or not self.old_path:
            return
        if sys.path != self.s.py_path + self.old_path:
            W = ("Could not deactivate python path for " 
                 "class CodeGenerator object")
            warnings.warn(W, UserWarning)
        else:
            del sys.path[:len(self.s.py_path)]

    def import_tables(self):
        finalize_parse_args(self.s)
        if getattr(self.s, "table_classes", None):
             return
        table_classes = self.s.table_classes = []
        for module, module_params in self.s.tables:
             m = import_module(module)
             table_class = m.Tables
             table_classes.append((table_class, module_params))

    def load_table_generator(self, table_generator, params):
        assert isinstance(params, ImmutableDict)
        assert isinstance(table_generator, TableGenerator)
        #tables = {}
        #directives = {}
        tables = self.s.table_classes
        mockup = self.s.mockup
        load_tables(table_generator, tables, params, mockup)

        """
        for table_class, table_params in self.s.table_classes:
             common_params = params.restrict(table_params)
             m_tables = table_class(**common_params)
             if mockup and getattr(m_tables, 'mockup_tables', None):
                 new_tables = m_tables.mockup_tables
             else:
                 new_tables = m_tables.tables
             if mockup and getattr(m_tables, 'mockup_directives', None):
                 new_directives = m_tables.mockup_directives
             else:
                 new_directives = m_tables.directives
             tables.update(new_tables)
             directives.update(new_directives)
             table_generator.set_tables(tables, directives)
        """

    def generate(self):
        finalize_parse_args(self.s)
        end_header = StringOutputFile()
        s = self.s
        out_dir = s.out_dir
        out_header = open_for_write(out_dir, s.out_header)
        write_warning_generated(out_header)
        tg = TableGenerator()
        if s.export_kwd:
            tg.export_kwd = s.export_kwd
        self.load_table_generator(tg, EmptyImmutableDict)
        if s.verbose:
            print("Generating header %s" % s.out_header)
        for param, c_files in s.c_files.items():
            self.load_table_generator(tg, param)
            for src_name, dest_name in c_files:
                src_file_name = self.find_file(src_name)
                src = open(src_file_name, "rt")
                if dest_name:
                    dest = open_for_write(out_dir, dest_name)
                    if s.verbose:
                        print("Generating file %s" % dest_name)
                    write_warning_generated(dest)
                    out_header.write("// %%%%FROM %s\n" % dest_name)
                    tg.generate(src, dest, out_header)
                    out_header.write("\n")
                    dest.close()
                else:
                     h_head, h_split, h_tail = split_header(src)
                     tg.generate(h_head, None, out_header, 'h')
                     end_header.write(h_split)
                     tg.generate(h_tail, None, end_header, 'h')

        end_header.copy_to(out_header)
        out_header.close()

    def c_files(self):
        finalize_parse_args(self.s)
        result = []
        for _, c_files in self.s.c_files.items():
            for _, dest_name in c_files:
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
                for src, dest in files:
                    print(" ", src, ",", dest)
        tables = getattr(s, 'tables', None)
        if tables:
            print("tables:")
            for table, par in tables:
                print(" ", table, par)

        attribs = [
            'source_header', 'out_header', 'py_path',
            'out_dir', 'export_kwd', 'verbose'
        ]
        for attr in attribs:
            print(attr + ':', getattr(s, attr, '<undefined>'))

           


def example():
    SAMPLE_ARGS = """ 
    --tables tt
    --source-header generate_code.py
    --source-path foo bar .
    --param p q 
    --tables t1 t2  
    --out-dir blah
    --set p 3
    --subst code code{p}
    --sources generate_code.py
    --set p 5
    --sources generate_code.py
    --out-dir source
    --out-header hhh
    --subst code code_again{p}
    --set p 3
    --sources generate_code.py
    --source-header generate_code.py
    --py-path shit
    """
    parser = generate_code_parser()
    s = parser.parse_args(SAMPLE_ARGS.split())
    print(s)
    print("\n")
    cg = CodeGenerator(s)
    cg.activate_py_path()
    sys.path.append("wtf")
    cg.deactivate_py_path()
    print("")
    parser.print_help()


if __name__ == "__main__":
    example()

