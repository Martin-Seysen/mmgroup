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
from collections import UserDict

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
from mmgroup.generate_c.generate_functions import eval_codegen_with
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



class NoDirectives:
    """This a fixed object ndicating that no directves are present"""
    pass



class TableGenerator(object):
    """Automatically generate a .c file and a .h file from a *source*
 
    Basic idea: We copy a *source* file to the .c file. The source 
    contains directives to enter precomputed tables, code  snippets 
    or constants into the .c file, and prototypes into the .h file.

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

             ``name`` : ``table`` .

        ``name`` is a string that is a valid python identifier.
        ``table`` is the python object to be referred by
        that name whenever an argument is evaluated as a
        python expression.

    :param directives:

        Here ``directives`` is a dictionary of shape::

            ``function_name`` : ``function``  .

        The function ``function`` is executed when a
        directive of shape::

             // %%function_name

        occurs as a directive in the *source* file. The arguments 
        following the directive are evaluated and passed to the 
        function ``function``. That function should return a
        string. This string is merged into the output .c file.

        If ``directives`` is ``mmgroups.generate_c.NoDirectives``
        then also built-in directives are not executed.

    :param verbose: 

        optional, may be set ``True`` for verbose console output


    External modules importing this class should use methods
    ``set_tables``, ``generate``, and ``table_size`` only.
    """

    m_kwd =  re.compile(r"\s*//\s*\%\%(\w+)(\*|\b)(.*)?")
    class m_no_kwds:
        def match(*args, **kwds):
            return None
    block_kwds = set(["FOR", "IF", "WITH"])
    ILLEGAL_INSIDE_BLOCK = "%s directive is illegal inside a codegen block"
    is_terminal = True, # True: Error if cannot valuate %{..}  expression
                        # False: Let expression unchanged in that case

    def __init__(self, tables = {}, directives={}, verbose = False):
        """Creates a Code generator with tables and directives
 
        """
        self.verbose = verbose
        self.names = {}       # Updated version of self.names
        self.start_sync()     # Reset all possibly pending items
        self.set_tables(tables, directives)

    def set_tables(self, tables = {}, directives = {}, args = ()):
        self.sync() 
        self.tables = tables
        self.reset_names()    # initialize self.names from self.tables
        # Enter user-defined directives into dictionary self.directives
        self.directives = {}
        if directives != NoDirectives:
            for fname, f in directives.items():
                if not isinstance(f, UserDirective):
                    try:
                        f = UserDirective(f)
                    except:
                        print(("\nError: Could not register directive"
                               + " %s for code generator\n") % fname )
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
             "PYX":          self.pyx_cmd,
             "FOR":          self.block_for,
             "WITH":         self.block_with,
             "JOIN":         self.do_join,
             "IF":           self.block_if,
        } )
 
        self.m_kwd = self.__class__.m_no_kwds        
        if directives != NoDirectives:
            self.directives.update(builtin_directives)
            self.m_kwd = self.__class__.m_kwd
            self.adjust_names()

        self.args = args     # positional arguments passed to formatting
                             # function, not used for C files

    def start_sync(self):
        """Reset all directives to the state 'completed'"""
        self.exported_table_names = {}
        #  dictionary python_table_name : c_table_name
        # 'python_table_name' is a key of dictionary tables.
        # 'c_table_name' is the corresponding name of that table
        # or entry in C.
        self.use_table_pending = False
        self.join_pending = False
        self.export_pending = 0
        self.export_args = ""
        self.export_kwd = ""
        self.block_stack = [] # stack for loop keywords such as FOR, IF, ..

        self.current_line = ""
        self.table_size_ = 0
        self.gen_c = True
        self.gen_h = False


    def sync(self):
        """Complain if there are any uncompleted directives

        This function also sets all relevant attriibutes such as
        expected when enetering a new source file.
        """
        ERR_DIRECTIVE = r"Directive %s has not been processed properly"
        if self.export_pending:
            raise ValueError(ERR_DIRECTIVE % "%%EXPORT")
        if self.use_table_pending:
            DIR = "%%EXPORT_TABLE or %%USE_TABLE"
            raise ValueError(ERR_DIRECTIVE % DIR)
        if len(self.block_stack):
            DIR = "%%" + self.block_stack[0]
            raise ValueError(ERR_DIRECTIVE % DIR)
        if self.join_pending:
            raise ValueError(ERR_DIRECTIVE % "%%JOIN")
        
        self.C_table_name = None            # C name of table pending
        self.C_table_export = False         # export C name of table pending if set
        self.exported_table_names.clear() 
        self.export_args = ""
        #self.export_kwd = ""
        self.current_line = ""
        self.gen_c = True
        self.gen_h = False


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
        try:
            name, table, format_  = eval_codegen_table(self.tables, args)
        except:
            _ , table, format_  = eval_codegen_table(self.names, args)
            name = ""
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
        self.use_table_pending = 1 + self.export_pending
        self.export_pending = 0
        self.export_args = "" 
        return "", ""

    def export_table(self, args, *_): 
        """built-in function EXPORT_TABLE"""
        self.use_table_pending = 2
        self.export_pending = 0
        export_str = self.export_kwd + "\n" if self.export_kwd else ""
        export_str_h = "// %%EXPORT_TABLE " + args + "\n"
        return export_str, export_str_h + export_str


    def comment(self, args, *_):
        """built-in function COMMENT, deprecated!!"""
        warnings.warn(
            r"The %%COMMENT directive in the code generator is deprecated!",
            DeprecationWarning
        )
        return "", ""


    def export_(self, args, *_):
        """built-in function EXPORT"""
        self.export_pending = True
        args = args.strip() 
        export_str = self.export_kwd + "\n" if self.export_kwd else ""
        export_str_h = "// %%EXPORT " + args + "\n"
        return export_str, export_str_h + export_str

    def set_export(self, args, *_):
        """built-in directive SET_EXPORT, deprecated!!!!"""
        err = "The SET_EXPORT directive is no longer supported"
        raise ValueError(err)


    def set_export_kwd(self, args, *_):
        """built-in function EXPORT_KWD"""
        self.export_kwd = args.strip()
        return "", ""


    def gen(self, args, *_):
        """built-in function GEN"""
        if len(self.block_stack):
            raise TypeError(self.ILLEGAL_INSIDE_BLOCK % "GEN")
        arg = direct_single_arg(args, "GEN") 
        self.gen_h =  arg.find('h') >= 0
        self.gen_c =  arg.find('c') >= 0
        self.gen_c =  self.gen_c or not (self.gen_c or self.gen_h)
        return "", ""
        

    def pyx_cmd(self, args, *_):
        """built-in function PYX"""
        err = "The PYX directive is no longer supported"
        raise ValueError(err)


    def docstr(self, args, *_):
        """built-in function PY_DOCSTR"""
        err = "The PY_DOCSTR directive is no longer supported"
        raise ValueError(err)


    def parse_block(self, kwd, source):
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
                level = len(self.block_stack) + len(kwd_stack) 
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
        raise TypeError("Missing END %s directive in codegen" % kwd)
        
      

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


    def block_for(self, args, quiet, source):
        """built-in function FOR"""
        local_vars, arg_list = eval_codegen_for(self.names, args)
        c_out, source, end_comment = self.parse_block("FOR", source)
        if quiet:
            c_out, end_comment = [], []
        h_out = []
        save_names = SaveDictFor(self.names, local_vars)
        self.block_stack.append("FOR")
        if self.join_pending:
            self.do_for_with_join(arg_list, source, save_names, c_out)
            self.join_pending = False
        else:
            for entry in arg_list:
                save_names.set_local(entry) # enter entry into self.names
                for c_1, h_1 in self.iter_generate_lines(iter(source)):
                    c_out.append(c_1)
                    h_out.append(h_1)
        self.block_stack.pop()
        save_names.restore()   # restore self.names
        return "".join(c_out + end_comment) , "".join(h_out)

    def block_with(self, args, quiet, source):
        """built-in function WITH"""
        local_vars, values = eval_codegen_with(self.names, args)
        c_out, source, end_comment = self.parse_block("WITH", source)
        if quiet:
            c_out, end_comment = [], []
        h_out = []
        save_names = SaveDictFor(self.names, local_vars)
        self.block_stack.append("WITH")
        save_names.set_local(values) # enter entry into self.names
        for c_1, h_1 in self.iter_generate_lines(iter(source)):
            c_out.append(c_1)
            h_out.append(h_1)
        self.block_stack.pop()
        save_names.restore()   # restore self.names
        return "".join(c_out + end_comment) , "".join(h_out)


    def do_join(self, args, quiet, source):
        """built-in function JOIN"""
        self.join_pending =  eval_codegen_join(self.names, args)
        return "", ""



    def parse_block_if(self, args, source):
        """Parser for reading an IF block from the source

        This is a variant of method parse_block() that also processes
        the 'ELSE' statementd appropriately.

        The function returns a triple

            (prefix_lines, code_lines, postfix_lines)

        similar to method parse_block(). List 'code_lines' contains
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
                level = len(self.block_stack) + len(kwd_stack) 
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


    def block_if(self, args, quiet, source):
        """built-in function IF"""
        c_out, source, end_comment = self.parse_block_if(args, source)
        if quiet:
            c_out, end_comment = [], []
        h_out = []
        self.block_stack.append("IF")
        for c_1, h_1 in self.iter_generate_lines(iter(source)):
            c_out.append(c_1)
            h_out.append(h_1)
        self.block_stack.pop()
        return "".join(c_out + end_comment) , "".join(h_out)


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
                    prototype = line[:found].strip()
                    h_out = "extern " + prototype +  ";\n"
                    self.C_table_export = True
        self.use_table_pending = 0
        return "", h_out

        

    def process_export(self, line):                
        """Do pending work of built-in function EXPORT"""
        if not self.export_pending:
            return
        prototype = line[:-1].strip()
        h_out = prototype + ";\n"
        self.export_pending = False
        return "", h_out

              
                   

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
        for l in source:
            if self.verbose:
                print("%d: %s" % (len(self.block_stack), l), end = "")
            m = self.m_kwd.match(l)
            self.current_line = l
            if m:
                fname, quiet, args = m.groups()
                if not quiet:
                    yield l, ""
                args = format_line(args, self.names, self.args) 
                c_out, h_out = self.directives[fname](args, quiet, source)
                self.adjust_names()      # a directive might mess up names
                yield c_out, h_out
            else:
                l = format_line(l, self.names, self.args,
                                           self.is_terminal)
                l = indent_subsequent_lines(l) 
                c_out = l if self.gen_c else "" 
                h_out = l if self.gen_h else ""
                yield c_out, h_out
                if self.export_pending:
                    yield self.process_export(l)
                if self.use_table_pending:
                    yield self.process_use_table(l)
                self.join_pending = False


    def iter_generate(self, source, gen='c'):
        """Iterator for generating .c and .h file form a source file.

        The function yield pairs of lines. The first line of each pair
        is to be written into the .c file and the second line is to
        be written into the .h file.

        Parameter source should be an iterator that yields the lines of
        the source file.
        """
        self.sync()
        gen = gen.strip()
        if gen:
            self.gen(gen)
        for c_out, h_out in self.iter_generate_lines(source):
            yield c_out, h_out
        self.sync()


    def generate(self, source, c_stream, h_stream = None, gen='c'):
        """Generate a .c and  a .h file form a source.

        Here ``source`` must be an iterator that yields the lines
        of the source file. This may be a readble object of class
        ``_io.TextIOWrapper`` as returned by the built-in function
        ``open()``.

        Parameters ``c_stream`` and ``h_stream`` must be either
        ``None`` or an instance of a class with a ``write`` method,
        e.g. a writable object of class ``_io.TextIOWrapper``.
        Then the output to the c file and to the h file is written
        to ``c_stream`` and to ``h_stream``. 
        """
        for c_out, h_out in self.iter_generate(source, gen):
            if c_out and c_stream:
                 c_stream.write(c_out) 
            if h_out and h_stream:
                 h_stream.write(h_out) 


    def table_size(self):
        return self.table_size_


    def display_export_file(self, text = None):
        """Display current export file at stdout"""
        if text: print(text)
        for line in self.export_file:
            print(" " + line)





def c_snippet(source, *args, **kwds):
    r"""Return a C code snippet as a string from a *source* string

    Here ``source`` is a string that is interpreted in the same way 
    as the text in a source file in method ``generate()`` of class 
    ``TableGenerator``.
    The function applies the code generator to the string ``source``
    and returns the generated C code as a string.

    All subsequent keyword arguments are treated as a dictionary and
    they are passed to the code generator in class ``TableGenerator``
    in the same way as parameter ``tables`` in the constructor of 
    that class. 

    One can also pass positional arguments to this function. In the
    ``source`` string they an be accessed as ``%{0}``, ``%{1}``, etc.

    A line starting with ``// %%`` is interpreted as a directive as
    in class ``TableGenerator``.

    The keyword ``directives`` is reserved for passing a dictionary
    of directives as in class ``TableGenerator``.

    In order to achieve a similar effect as generating C code with::

         tg = TableGenerator(tables, directives)
         tg.generate(source_file, c_file)

    you may code e.g.::
        
        c_string = c_snippet(source, directives=directives, **tables) 

    If this function cannot evaluate an expression of shape ``%{xxx}`` 
    then the expression it is not changed; so a subsequent code 
    generation step may evaluate that expression. An unevaluated 
    argument in a directive leads to an error.
    """ 
    src = source.splitlines(True)
    try:
        directives = kwds['directives']
    except KeyError:
        directives = {}
    tg = TableGenerator()
    tg.is_terminal = False
    tg.set_tables(kwds, directives, args)
    return "".join(c for c, h in tg.iter_generate(src)) 



def make_doc():
    raise NotImplementedError


class TableGeneratorStream(TableGenerator):
    """This is a deprecated version of class TableGenerator

    It is used for compatibility with the old code generation process
    before switching to meson.

    It will contain the methods of class TableGenerator that
    will be dprecated when bulding the project with meson.
    """
    def __init__(self, *args, **kwds):
        super(TableGeneratorStream, self).__init__(*args, **kwds)

    def out_auto_headlines(self):
        s = """{1}
// This {0} file has been created automatically. Do not edit!!!
{1}

"""
        cmt = "/" * 77
        return s.format("C", cmt), s.format("C header", cmt)



    def iter_source(self, source):
        """An iterator that yields the lines of the source file.

        The source file (or the source string) is given by parameter 
        *source*. Description of parameter *source* see method 
        generate().

        The lines yielded by this function are fed into method
        iter_generate() of this class.
        """
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
               Then that string is interpreted as a name of a file

             * A list of strings.
               This is interpreted as a sequence of strings or file names
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
            c_file.write(c_out)
        if  not h_filename or isinstance(h_filename, file): 
            h_file = h_filename                # output .h file
        else:
            h_file = open(h_filename,"wt")     # output .h file
            h_file.write(h_out)

        super(TableGeneratorStream, self).generate(inp, c_file, h_file)

        if c_filename and not isinstance(c_filename, file):
            c_file.close()
        if h_filename and not isinstance(h_filename, file):
            h_file.close()

