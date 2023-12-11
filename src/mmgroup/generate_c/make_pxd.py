"""Auxiliary functions for writing tables into C program
"""

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



import sys
import re
import os
import warnings
import ast
from io import IOBase


#from mmgroup.generate_c.make_c_tables import TableGenerator

m_kwd =  re.compile(r"\s*//\s*\%\%(\w+)\s+(.*)")
m_function = re.compile(r"[^()]+\([^()]*\)")
m_voidargs_function = re.compile(r"[^()]+\(\s*void\s*\)")




def iter_exports_from_header(file):
    """Yield list of export tuples from a generated header file

    Here ``file`` is the name of the heder file to be read.

    The function yields triples ``(directive, args, prototype)``
    describing function to be exported, e.g. via a .pxd file.
    Here ``directive`` is the the directive read from the header
    file. At present ``directive`` will be one of the strings
    'EXPORT',  'EXPORT_TABLE', or 'FROM'. Component ``args``
    is a single string containing the arguments read from the
    directive. Component ``prototype`` is the prototype for the
    exported function or table, such that it can be entered into
    a .pxd file to be generated. For a 'FROM' directive, component
    ``args`` is the name of the C source file from which the
    subsequent functions and tables are exported; and ``prototype``
    is the empty string.
    """
    f = open(file, "rt")
    export_pending = None
    def prototype(line):
        if line.endswith(';'):
            line = line[:-1].rstrip()
        if m_voidargs_function.match(line):
            line = line[:line.index('(')] + "()"
        return line

    for line in f:
        line = line.strip()
        m = m_kwd.match(line)
        if m:
            kwd, args = [s.strip() for s in m.groups()] 
            if kwd in ['EXPORT', 'EXPORT_TABLE']:
                export_pending = kwd, args
            elif kwd == 'FROM':
                yield kwd, args, ""
        elif export_pending:
            kwd, args = export_pending
            if kwd == 'EXPORT': 
                 if m_function.match(line):
                     yield kwd, args, prototype(line)
                     export_pending = None
            elif kwd == 'EXPORT_TABLE': 
                 if line.startswith('extern'):
                     yield kwd, args, prototype(line)
                     export_pending = None
            else:
                 export_pending = None
            if line.isspace():
                 export_pending = None




def iter_pxd_in(source):
    """An iterator that yields the lines of the .pxd file
    of a source parameter.

    The source file (or the source string) is given by parameter 
    *pxd_in*. Description of parameter *pxd_in* in function
    generate_pxd().
    """
    if isinstance(source, IOBase):
        for l in source:
            yield l
    else:
        try:
            is_filename = not "\n" in source
        except:
            is_filename = True
        if is_filename:
            with open(source, "rt") as f:
               for l in f:
                   yield l
        else:
            if source[-1] == "\n":
                 source = source[:-1]
            for l in source.split("\n"):
                 yield l + "\n"     



def make_pxd_entry(kwd, par, prototype):
    s = ""
    if 'p' in par:
        if kwd == 'EXPORT' and 'x' in par:
            s += "    # PYX <wrap pxi>\n"
        if kwd in ['EXPORT', 'EXPORT_TABLE']:
            s += "    %s\n" % prototype
        if kwd == 'FROM':
            s = "     # from " + par + "\n"
    return s  
    





def pxd_from_h(pxd_out, h_in, pxd_in = None, h_name = None, nogil = False):
    r"""Create a ``.pxd`` file from a C header file.

    Here parameter ``h_in`` is the name of a C header file that has
    usually been generated with  method ``generate`` of
    class ``TableGenerator``. In such a header file, some exported
    function are preceded by a comment of shape

        ``// %%EXPORT <parameter>``

    There ``<parameter>`` is a string of letters describing the way how
    a function is exported. The function is entered into the output 
    ``.pxd`` file if ``<parameter>`` contains the letter ``'p'``.

    Method ``generate_pxd`` creates a .pxd file with name ``pxd_file`` 
    that contains information about these exported function.

    A .pxd file is used in Cython. Its content is::

       cdef extern from <h_file_name>:
           <exported function 1>
           <exported function 2>
           ...

    Here ``<h_file_name>`` is given by parameter ``h_name``.

    :param pxd_out:

        Name of the output .pxd file.

    :param h_in:

        Name of the header file from which external functions are copied
        into the .pxd file as a ``cdef extern from`` statement.

    :param pxd_in:
 
         An optional name of a input .pxd file to be copied in the
         output .pxd file in front of the ``cdef`` statement. If
         that parmeter is a string containing a newline character
         ``\n`` then that string is copied directly into the
         output .pxd file.

    :param h_name:
         
         Name of .h file to to be placed into the statement
         ``cdef extern from <h_file_name>``.

    :param nogil:

         Optional, exported functions are declared as ``nogil``
         when set.
    """
    s_gil = " nogil" if nogil else ""
    must_open_pxd_out = not isinstance(pxd_out, IOBase)
    if must_open_pxd_out: 
        f = open(pxd_out, "wt")
    else:
        f = pxd_out
    print("# This .pxd file has been generated automatically. Do not edit!\n",
        file = f)
    if pxd_in:
        for l in iter_pxd_in(pxd_in):
           print(l, file = f, end = "")
        print("", file = f)
    if h_name is None:
        _, h_name = os.path.split(h_in)
    print('cdef extern from "{0}"{1}:'.format(h_name, s_gil), file = f)
    for directive, args, prototype in iter_exports_from_header(h_in):
        s = make_pxd_entry(directive, args, prototype)
        print(s, file = f, end = "")
    if must_open_pxd_out: 
        f.close()

generate_pxd = pxd_from_h



def c_file_list_from_h(h_in):
    c_list = []
    for directive, args, _ in iter_exports_from_header(h_in):
        if directive == 'FROM':
            c_list.append(args)
    return c_list


if __name__ == "__main__":
    h_in = os.path.join('..', 'dev', 'c_files', 'clifford12.h')
    generate_pxd(sys.stdout, h_in, nogil=True)





