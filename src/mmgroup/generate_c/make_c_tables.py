"""Auxiliary functions for writing tables into C program"""

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import types
import sys
import re
import os
import warnings
import itertools
import ast
from io import IOBase
from operator import __or__, __xor__, __and__
from numbers import Integral
from collections.abc import Iterable
from functools import partial
import numpy






if sys.version_info[0] >= 3:
    file = IOBase 
    long = int    


import mmgroup.generate_c.generate_functions 


from mmgroup.generate_c.generate_functions import isinteger
from mmgroup.generate_c.generate_functions import format_item
from mmgroup.generate_c.generate_functions import sizeof_format
from mmgroup.generate_c.generate_functions import make_table
from mmgroup.generate_c.generate_functions import py_doc_to_comment
from mmgroup.generate_c.generate_functions import format_line
from mmgroup.generate_c.generate_functions import line_as_comment
from mmgroup.generate_c.generate_functions import named_table
from mmgroup.generate_c.generate_functions import eval_codegen_args
from mmgroup.generate_c.generate_functions import eval_codegen_table
from mmgroup.generate_c.generate_functions import eval_codegen_for
from mmgroup.generate_c.generate_functions import eval_codegen_else
from mmgroup.generate_c.generate_functions import SaveDictFor
from mmgroup.generate_c.generate_functions import safe_eval
from mmgroup.generate_c.generate_functions import direct_single_arg
from mmgroup.generate_c.generate_functions import prepend_blanks
from mmgroup.generate_c.generate_functions import indent_subsequent_lines
from mmgroup.generate_c.generate_functions import eval_codegen_join
from mmgroup.generate_c.generate_functions import UserDirective
from mmgroup.generate_c.generate_functions import UserFormat
from mmgroup.generate_c.generate_functions import built_in_formats

class TableGenerator(object):
    """Automatically generate a .c file and a .h file from a *source*
 
    Basic idea: We copy a *source* file to the .c file. The source 
    contains directives to enter precomputed tables, code  snippets 
    or constants into the .c file, and protoypes into the .h file.

    We copy the content of the *source* file directly to the .c
    file.   In the *source* file we also parse certain directives 
    and we replace them by automatically generated code snippets 
    in the .c file and in the .h file. 

    The arguments given to the constructor of this class specify
    the meaning of the directives.

    For a more detailed description of these directives see module
    ``make_c_tables_doc.py`` 

    :param tables:

        Here ``tables`` is a dictionary of with entries of shape:: 

             name : table .

        ``name`` is a string that is a valid python identifier.
        ``table`` is the python object to be referred by
        that name whenever an argument is evaluated as a
        python expression.

    :param directives:

        Here ``directives`` is a dictionary of shape::

            function_name : function  .

        The function ``function`` a directive is executed when a
        direcrive of shape::

             // %%function_name

        occurs as keyword in a directive in the *soure* file.
        The following arguments are evaluated and passed to  
        function ``function``. The function should return a
        string. This string is merged into the output .c file.

    :param verbose: 

        optional, may be set ``True`` for verbose console output
    """

    m_kwd =  re.compile(r"\s*//\s*\%\%(\w+)(\*|\b)(.*)?")
    block_kwds = set(["FOR", "IF"])
    ILLEGAL_INSIDE_BLOCK = "%s directive is illegal inside a codegen block"


    def __init__(self, tables, directives={}, verbose = False):
        """Creates a Code generator with tables and directives
 
        """
        self.verbose = verbose
        self.tables = tables  # This is never changed
        self.names = {}       # Updated version of self.names
        self.reset_names()    # initialize self.names from self.tables
        self.exported_table_names = {}
        #  dictionary python_table_name : c_table_name
        # 'python_table_name' is a key of dictionary tables.
        # 'c_table_name' is the corresponding name of that table
        # or entry in C.

        # Enter user-defined directives into dictionary self.directives
        self.directives = {}
        for fname, f in directives.items():
            if not isinstance(f, UserDirective):
                try:
                    f = UserDirective(f)
                except:
                    print("\nError:")
                    print("Could not register directive %s for code generator\n" % fname)
                    raise
            f.set_code_generator(self)
            self.directives[fname] = f
         

        # Enter built-in directives into dictionary self.directives
        # Dictionary 'built_in_directives' of is shape
        #     name : <directive>.
        # Here <directives> is a member functions of this class 
        builtin_directives = ( {
             "TABLE":        self.table, 
             "USE_TABLE":    self.use_table,
             "EXPORT":       self.export_,
             "SET_EXPORT":   self.set_export,
             "EXPORT_KWD":   self.set_export_kwd,
             "EXPORT_TABLE": self.export_table,
             "GEN":          self.gen,
             "COMMENT":      self.comment,
             "PY_DOCSTR":    self.docstr,
             "PYX":          self.pyx,
             "FOR":          self.loop_for,
             "JOIN":         self.do_join,
             "IF":           self.loop_if,
        } )
        self.directives.update(builtin_directives)

        self.use_table_pending = False
        self.export_pending = 0
        self.pxd_export_pending = 0
        self.export_kwd = ""
        self.args = ()       # positional arguments passed to formatting
                             # function, not used for C files
        self.loop_stack = [] # stack for loop keywords such as FOR, IF, ..

        self.source_name = None  # name of C skeleteon file
        self.current_line = ""


    def reset_names(self):
        """Reset dictionary self.names to its original value self.tables

        This means that all C names for tables generated by the
        %%TABLE directive will be lost.
        """
        self.names.clear()
        #names.update(safe_locals)
        self.names.update(built_in_formats)
        self.names.update(self.tables)
        self.names.update({
            "NAMES" : self.names,  
            "TABLES" : self.tables, 
        })
        self.adjust_names()

    def adjust_names(self):
        for _, entry in self.names.items():
            if isinstance(entry, UserDirective):
                 entry.set_code_generator(self)


    def _check_names(self):
        for _, entry in self.names.items():
            if isinstance(entry, UserDirective):
                 assert entry.tg == self


    def table(self, args, *_):
        """built-in function TABLE"""
        name, table, format_  = eval_codegen_table(self.tables, args)
        if name and self.C_table_name:
            self.names[name] = named_table(self.C_table_name, table)
            if self.C_table_export:
                self.exported_table_names[name] = self.C_table_name
        self.C_table_name = None
        self.C_table_export = False
        try:
            self.table_size_ += len(table) * sizeof_format(t_args[0])
        except:
            pass
        return make_table(table, format_), ""

    def use_table(self, args, *_):
        """built-in function USE_TABLE"""
        if len(self.loop_stack):
            raise TypeError(self.ILLEGAL_INSIDE_BLOCK % "USE_TABLE")
        self.use_table_pending = 1 + self.export_pending
        self.export_pending = 0
        self.pxd_export_pending = 0
        return "", ""

    def export_table(self, args, *_): 
        """built-in function EXPORT_TABLE"""
        if len(self.loop_stack):
            raise TypeError(self.ILLEGAL_INSIDE_BLOCK % "EXPORT_TABLE")
        self.use_table_pending = 2
        self.export_pending = 0
        self.pxd_export_pending = 0
        return self.export_kwd, self.export_kwd


    def comment(self, args, *_):
        """built-in function COMMENT"""
        return "", ""


    def export_(self, args, *_):
        """built-in function EXPORT"""
        if len(self.loop_stack):
            raise TypeError(self.ILLEGAL_INSIDE_BLOCK % "EXPORT")
        self.export_pending = True
        self.pxd_export_pending = len(args) and 'p' in args
        return self.export_kwd, self.export_kwd

    def set_export(self, args, *_):
        """built-in function SET_EXPORT, deprecated!!!!"""
        raise ValueError("Function no longer supported")
        if len(self.loop_stack):
            raise TypeError(self.ILLEGAL_INSIDE_BLOCK % "SET_EXPORT")
        #args = format_line(args, self.names, self.args)
        args = tuple(s.strip() for s in args.split("=")) 
        name, value = args[0], args[1]
        self.exported_table_names[name] = value
        self.names[name] = value
        return "", ""


    def set_export_kwd(self, args, *_):
        """built-in function EXPORT_KWD"""
        self.export_kwd = args.strip() + "\n"
        return "", ""


    def gen(self, args, *_):
        """built-in function GEN"""
        if len(self.loop_stack):
            raise TypeError(self.ILLEGAL_INSIDE_BLOCK % "GEN")
        arg = direct_single_arg(args, "GEN") 
        self.gen_h =  arg.find('h') >= 0
        self.gen_c =  arg.find('c') >= 0
        self.gen_c =  self.gen_c or not (self.gen_c or self.gen_h)
        return "", ""
        

    def pyx(self, args, *_):
        """built-in function PYX"""
        if len(self.loop_stack):
            raise TypeError(self.ILLEGAL_INSIDE_BLOCK % "PYX")
        #args = format_line(args, self.names, self.args)
        self.pxd_entries.append("# PYX " + args.strip())
        return "", ""


    def docstr(self, args, *_):
        """built-in function PY_DOCSTR"""
        args= eval_codegen_args(self.tables, args)
        entry = args[0]
        if not isinstance(entry, str):
            entry = entry.__doc__
        if len(args) > 1 and args[1]:
            entry = format_line(entry, self.names, self.args)
        return py_doc_to_comment(entry), ""



    def parse_loop(self, kwd, source):
        """General parser for reading a block from the source 

        Block must nested in the form:

        // %%KWD1
           ...
           // %%KWD2
           ...
           // %%END KWD2
           ...
        // %%END KWD1

        with indentation ignored. Keywords KWD1, KWD2 etc. are block
        keywords. The set of  block keywords is given by
        self.block_kwds.

        Once a directive with a block keyword has been parsed, a parser
        for reading the complete block corresponding to that block 
        keyword should be called. So in the example above, this parser
        reads the lines up to and including the final line
        '// %%END KWD1' of te block. Note that the first line
        '// %%KWD1' of the block has already ben read by the function
        calling this method. 

        'kwd' is the block keyword already read by the calling funtion.
        *source* is the iteratornthat provides the input lines.

        The function returns a triple

            (prefix_lines, code_lines, postfix_lines)

        of list of strings representing input lines. The lines in the
        list 'code_lines' form the body of the block, which will be
        processed for code generation by a subsequent function.

        'prefix_lines' and 'postfix_lines' are lists of comment lines 
        to be copied from the soure to the generated code file.

        Here in the standard case, 'prefix_lines' is empty and 
        'postfix_lines' contains the END statement of the block.
        """
        kwd_stack = [kwd]
        out_source = []
        for l in source:
            if self.verbose:
                level = len(self.loop_stack) + len(kwd_stack) 
                print("%d. %s" % (level, l), end = "")
            out_source.append(l)
            m = self.m_kwd.match(l)
            if m:
                sub_kwd, _, args = m.groups()
                if sub_kwd in self.block_kwds:
                    kwd_stack.append(sub_kwd)
                elif sub_kwd == "END":
                    end_kwd = direct_single_arg(args, kwd = "END")
                    match_kwd = kwd_stack.pop()
                    if end_kwd != match_kwd:
                        kw = ("END " + end_kwd, match_kwd)
                        s =" % does not match % in codegen directive"
                        raise TypeError(s % kw) 
                    if len(kwd_stack) == 0:
                        return [], out_source[:-1], out_source[-1:]
        raise TypeError("Missing END %d directive in codegen" % kwd)
        
      

    def do_for_with_join(self, arg_list, source, save_names, c_out):
        """Special inner loop of FOR directive with JOIN prefix"""
        infix, suffix = self.join_pending
        for entry in arg_list:
            save_names.set_local(entry) # enter entry into self.names
            for c_out1, _ in self.iter_generate_lines(iter(source)):
                c_out.append(c_out1)
            if len(c_out) and  c_out[-1][-1:] == "\n":
                c_out[-1]  = c_out[-1][:-1] 
            c_out.append(infix)
        if len(c_out):
            c_out[-1] = suffix


    def loop_for(self, args, quiet, source):
        """built-in function FOR"""
        local_vars, arg_list = eval_codegen_for(self.names, args)
        c_out, source, end_comment = self.parse_loop("FOR", source)
        if quiet:
             c_out, end_comment = [], []
        save_names = SaveDictFor(self.names, local_vars)
        self.loop_stack.append("FOR")
        if self.join_pending:
            self.do_for_with_join(arg_list, source, save_names, c_out)
            self.join_pending = False
        else:
            for entry in arg_list:
                save_names.set_local(entry) # enter entry into self.names
                for c_out1, _ in self.iter_generate_lines(iter(source)):
                    c_out.append(c_out1)
        self.loop_stack.pop()
        save_names.restore()   # restore self.names
        return "".join(c_out + end_comment) , ""


    def do_join(self, args, quiet, source):
        """built-in function JOIN"""
        self.join_pending =  eval_codegen_join(self.names, args)
        return "", ""



    def parse_loop_if(self, args, source):
        """Parser for reading an IF block from the source

        This is a variant of method parse_loop() that also processes
        the 'ELSE' statementd appropriately.

        The function returns a triple

            (prefix_lines, code_lines, postfix_lines)

        similar to method parse_loop(). List 'code_lines' contains
        only the lines not excluded by IF or ELSE clauses. List
        'prefix_lines' contains the ELSE clause after which valid
        lines follow (if any).
        """ 
        kwd_stack = ["IF"]
        out_prefix = []
        out_source = []
        else_ok = True
        else_done = condition = bool(safe_eval(args, self.names))
        for l in source:
            if self.verbose:
                level = len(self.loop_stack) + len(kwd_stack) 
                print("%d. %s" % (level, l), end = "")
            new_condition = condition
            m = self.m_kwd.match(l)
            if m:
                sub_kwd, _, args = m.groups()
                if sub_kwd in self.block_kwds:
                    kwd_stack.append(sub_kwd)
                elif sub_kwd == "END":
                    end_kwd = direct_single_arg(args, kwd = "END")
                    match_kwd = kwd_stack.pop()
                    if end_kwd != match_kwd:
                        kw = ("END " + end_kwd, match_kwd)
                        s =" % does not match % in codegen directive"
                        raise TypeError(s % kw) 
                    if len(kwd_stack) == 0:
                        if not else_done:
                            out_prefix = []
                        return out_prefix, out_source, [l]
                elif sub_kwd == "ELSE" and len(kwd_stack) == 1:
                    new_condition = False
                    if not else_ok:
                        raise TypeError("Misplaced ELSE in IF directive in codegen")
                    if not else_done:
                        out_prefix = [l]
                        args = args.format(*self.args, **self.names) 
                        new_condition, else_ok = eval_codegen_else(
                            self.names, args)
                        else_done = new_condition
            if condition and new_condition:
                out_source.append(l)
            condition = new_condition 
             
        raise TypeError("Missing END IF directive in codegen")


    def loop_if(self, args, quiet, source):
        """built-in function IF"""
        c_out, source, end_comment = self.parse_loop_if(args, source)
        if quiet:
            c_out, end_comment = [], []
        self.loop_stack.append("IF")
        for c_out1, _ in self.iter_generate_lines(iter(source)):
             c_out.append(c_out1)
        self.loop_stack.pop()
        return "".join(c_out + end_comment) , ""


    def process_use_table(self, line):
        """Do pending work of built-in function USE_TABLE"""
        if not self.use_table_pending:
            return
        self.C_table_name, h_out = False, ""
        found = line.find("[")
        if found >= 0:
            l1 = line[:found].split(" ")
            self.C_table_name = l1[-1].strip()
            if self.use_table_pending > 1:
                found = line.find("=")
                if found:
                    h_out = "extern " + line[:found].strip() + ";\n"
                    self.C_table_export = True
        self.use_table_pending = 0
        return "", h_out

        

    def process_export(self, line):                
        """Do pending work of built-in function EXPORT"""
        if not self.export_pending:
            return
        h_out =  line[:-1] + ";\n"
        if self.pxd_export_pending:
            self.pxd_entries.append(line[:-1].strip())
        self.export_pending = False
        self.pxd_export_pending = False
        return "", h_out

              
       

    def out_auto_headlines(self):
        s = """{1}
// This {0} file has been created automatically. Do not edit!!!
{1}

"""
        cmt = "/" * 77
        return s.format("C", cmt), s.format("C header", cmt)

            

    def iter_generate_lines(self, source):
        """This iterator is the workhorse for method generate(line).

        It reads input from iterator source, and it yields pairs
        (c_line, h_line), where c_line is written to the generated
        .c file and h_line is written to the generated .h file.

        This method may is also be called by block processing directives
        corresponding e.g. to an IF or FOR directive. In case of a
        FOR directive, this method called once for each iteration 
        of the FOR loop.
        """
        self.adjust_names()
        self.export_pending = False         # set when EXPORT action pending
        self.pxd_export_pending = False     # set when EXPORT include .pxd file
        self.use_table_pending = False      # set when USE_TABLE action pending
                                            # 1 = just remember C name of table
                                            # 2 = also create prototype
        self.join_pending = False           # set when JOIN action is pending
        for l in source:
            if self.verbose:
                print("%d: %s" % (len(self.loop_stack), l), end = "")
            m = self.m_kwd.match(l)
            self.current_line = l
            if m:
                fname, quiet, args = m.groups()
                if not quiet:
                    yield l, ""
                args = args.format(*self.args, **self.names)
                c_out, h_out = self.directives[fname](args, quiet, source)
                self.adjust_names()        # a directive might mess up names
                yield c_out, h_out
            else:
                #self._check_names()
                l = format_line(l, self.names, self.args)
                l = indent_subsequent_lines(l) 
                c_out = l if self.gen_c else line_as_comment(l)
                h_out = l if self.gen_h else ""
                yield c_out, h_out
                if self.export_pending:
                    yield self.process_export(l)
                if self.use_table_pending:
                    yield self.process_use_table(l)
                self.join_pending = False
        #self.reset_names()


    def iter_generate(self, source):
        """Iterator for generating .c and .h file form a source file.

        The function yield pairs of lines. The first line of each pair
        is to be written into the .c file and the second line is to
        be written into the .h file.

        Parameter source should be an iterator that yields the lines of
        the source file.
        """
        self.C_table_name = None            # C name of table pending
        self.C_table_export = False         # export C name of table pending if set
        self.gen_c = True                   
        self.gen_h = False
        self.exported_table_names.clear()
        self.table_size_ = 0
        self.pxd_entries = []

        for c_out, h_out in self.iter_generate_lines(source):
            yield c_out, h_out



 
    @staticmethod
    def print_file(text, file):
        if text and file:
            print(text, end = "", file = file) 



    def iter_source(self, source):
        """An iterator that yields the lines of the source file.

        The source file (or the source string) is given by parameter 
        *source*. Description of parameter *source* see method 
        generate().

        The lines yielded by this function are fed into method
        iter_generate() of this class.
        """
        self.source_name = None
        if isinstance(source, list):
            for i, f in enumerate(source):
                for l in self.iter_source(f):
                    yield(l)
                if i < len(source) - 1:
                    yield "\n"
        elif isinstance(source, file):
            self.gen_c = True                   
            self.gen_h = False
            for l in source:
                yield l
        else:
            self.gen_c = True                   
            self.gen_h = False
            try:
                is_filename = not "\n" in source
            except:
                is_filename = True
            if is_filename:
                if not self.source_name:
                    self.source_name = source
                if self.verbose:
                    print(" > Opening file %s" % source)
                with open(source, "rt") as f:
                    for l in f:
                        yield l
                if self.verbose:
                    print(" > Closing file %s" % source)
               
            else:
                if source[-1] == "\n":
                    source = source[:-1]
                for l in source.split("\n"):
                    yield l + "\n"     
                     

    def generate(self, source, c_filename, h_filename=None):
        r"""Generate .c file and .h file from source file

        :param source:

           The source file to processed for code generation.
           This may be

             * A string containing at least one *newline* character.
               Then that string is directly interpreted as the content
               of the source.

             * Any other string. 
               Then that string is intepreted as a name of a file

             * A list of strings.
               This is intepreted as a sequence of strings or file names
               according to the last two rules. Then the items given by
               the entries of that sequence are concatenated.      

        :param c_filename:

           Name of the .c file to be generated or ``None``. 
           The extension must be part of that name.
      
        :param h_filename:

           Name of the .h file to be generated or ``None``
           (default). 
           The extension must be part of that name.
      
        """
        inp = self.iter_source(source)         # input source file
        c_out, h_out = self.out_auto_headlines()
        if not c_filename or isinstance(c_filename, file):
            c_file = c_filename                # output .c file
        else:
            c_file = open(c_filename,"wt")     # output .c file
            self.print_file(c_out, c_file)
        if  not h_filename or isinstance(h_filename, file): 
            h_file = h_filename                # output .h file
        else:
            h_file = open(h_filename,"wt")     # output .h file
            self.print_file(h_out, h_file)

        for c_out, h_out in self.iter_generate(inp):
            self.print_file(c_out, c_file)
            self.print_file(h_out, h_file)
        inp = None

        if c_filename and not isinstance(c_filename, file):
            c_file.close()
        if h_filename and not isinstance(h_filename, file):
            h_file.close()

    def table_size(self):
        return self.table_size_


    def generate_pxd(self, pxd_file, h_file_name, source = None, nogil = False):
        """Create .pxd file 'pxd_file' from last call to generate().

        After calling method generate(), all exported functions 
        preceeded by an ``EXPORT`` directive are collected automatically. 

        Method ``generate_pxd`` creates a .pxd file with name ``pxd_file`` 
        that contains information about these exported function.

        A .pxd file is used in Cython. Its content is::

           cdef extern from <h_file_name>:
               <exported function 1>
               <exported function 2>
               ...

        The name  of the .h file is given as parameter and the exported
        functions listed below that statement are those collected by
        the last recent call to method generate(). 
 
        :param pxd_file:

            Name of the output file.

        :param h_file_name:

            Name of the .h file to be placed in the
            ``cdef extern from`` statment.
    
        :param source:
     
             An optional string written into the output file
             in front of the ``cdef`` statement.  Format of parameter 
             *source* is the same as in method generate()

        :param nogil:

             Optional, exported functions are declared as ``nogil``
             when set.
        """
        s_gil = " nogil" if nogil else ""
        f = open(pxd_file, "wt")
        print("# This .pxd file has been generated automatically. Do not edit!\n",
            file = f)     
        if source:
            old_source_name = self.source_name
            for l in self.iter_source(source):
                print(l, file = f, end = "")
            print("", file = f)
            self.source_name = old_source_name
        print('cdef extern from "{0}"{1}:'.format(h_file_name, s_gil), file = f)  
        for l in self.pxd_entries:
            print("    " + l, file = f)
        if len(self.pxd_entries) == 0:
            print("    pass", file = f)
        f.close()








def make_doc(source_file, output_file, tables = None):
    r"""Extract documentation from C file generated by class ``TableGenerator``

    Here parameter ``source_file`` describes a source file to be used
    as input for an instance of class ``TableGenerator``. 

    Parameter ``source_file may`` be:
        * the name of a C source file
        * a string containing C statements with at least one *newline*
          character
        * a list of file names or string as defined above

    Parameter ``output_file`` is the name of an output file to be generated.
    The output file contains function declarations and comments for
    each in the same way as in the processed source file. The bodies
    of the C functions are dropped.

    The source file is processed essentially in the same way as a
    source file is procesed by an instance of class  ``TableGenerator``.

    If one of the ``TableGenerator`` directives ``EXPORT`` or ``COMMENT``
    is found in the source file then everything of the source file is 
    copied into the  output file up to and including the first comment.
 
    Here a comment line must begin with ``/*`` (preceeded by white
    space only) and end with ``*/``, or it must be as sequence of lines
    containing  ``//`` preceeded by white  space only. 

    Parameter ``tables`` is an optional dictionary that plays the same 
    role as in the constructor of class TableGenerator.
    """
    def iter_source(source_file):
        if isinstance(source_file, list):
            for i, f in enumerate(source_file):
                for l in self.iter_source(f):
                    yield(l)
                if i < len(source_file) - 1:
                    yield "\n"
        else:
            try:
                is_filename = not "\n" in source
            except:
                is_filename = True
            if is_filename:
                with open(source_file, "rt") as f:
                    for l in f:
                        yield l               
            else:
                if source_file[-1] == "\n":
                    source_file = source_file[:-1]
                for l in source_file.split("\n"):
                    yield l + "\n"     


    g = open(output_file, "wt")
    m_copy_comment =   re.compile(r"\s*//\s*\%\%(EXPORT|COMMENT)\b")
    m_begin_comment =  re.compile(r"\s*/\*.*(\*/)?")
    m_end_comment = re.compile(r".*\*/")
    m_is_comment = re.compile(r"\s*//")
    m_stop = re.compile(r"\s*[{$]")
    m_any_directive =   re.compile(r"\s*//\s*\%\%[A-Za-z0-9_]*\b")
    comment_pending = in_comment = copy = is_keyword = False 
    for l in iter_source(source_file):
        is_keyword = False
        if  m_copy_comment.match(l):
            comment_pending  = is_keyword = True
            in_comment, copy = False, True
        elif m_begin_comment.match(l):
            if comment_pending:
                copy, commment_pending = True, False
                closing = m_begin_comment.match(l).groups()[0]
                in_comment = not closing
                comment_pending = False
        elif in_comment and m_end_comment.match(l):
            in_comment = False
        elif m_is_comment.match(l):
            if comment_pending:
                copy, comment_pending = True, False
        else:
            if not in_comment:
                if copy and not comment_pending:
                    print("\n\n",  file = g)
                    copy = False

        if  copy and not is_keyword:
            try:
                l = format_line(l, tables)
            except:
                pass 
            if not m_any_directive.match(l): 
                print(l, end = "", file = g)
    g.close()        
 



def c_snippet(source, *args, **kwds):
    r"""Return a C code snippet as a string from a *source* string

    Here *source* is a string that is interpreted in the same way as
    a source file in method ``generate()`` of class TableGenerator.

    All subsequent keyword arguments are treated as a dictionary and
    they are passed to the code generator in class TableGenerator in 
    the same way as parameter 'tables' in the constructor of that
    class.

    Thus ``c_snippet(source, x, a = 'b')`` has pretty much the same 
    effect ``as source.format(x, a = 'b')``. This means that 
    positional and keyword arguments of this function are eventually 
    passed to to the standard string formatting method ``.format()``,
    with the following exceptions: 
 
       - A opening curly bracket not followed by a closing one or a 
         closing curly bracket not preceeded by an opening one is not 
         passed to the formatting method. So we may write the string  
         *source*  almost in usual C syntax, but expressions in curly 
         braces inside a line are passed to the formatting method.

       - A line starting with ``// %%`` is interpreted as in class
         TableGenerator.

       - keyword ``directives`` is reserved for passing a list of
         directives as in class TableGenerator.

    In order to achieve a similar effect as generating C code with::

         tg = TableGenerator(tables, directives)
         tg.generate(source_file, c_file)

    you may code e.g.::
        
        c_string = c_snippet(source, directives=directives, **tables) 
    """ 
    src = source.splitlines(True)
    try:
        directives = kwds['directives']
    except KeyError:
        directives = {}
    tg = TableGenerator(kwds, directives)
    tg.args = args
    return "".join(c for c,h in tg.iter_generate(src)) 



