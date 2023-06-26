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


from mmgroup.generate_c.make_c_tables import TableGenerator



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


def iter_export_in(source):
    """An iterator that yields the lines of the export file
    of a source parameter.

    The source file (or the source string) is given by parameter 
    *export_in*. Description of parameter *export_in* in function
    generate_pxd().
    """
    if isinstance(source, IOBase):
        for l in source:
            yield l.strip()
    elif isinstance(source, TableGenerator):
        #source.display_export_file("TableGenerator")
        for l in source.export_file:
            yield l.strip() 
    else:
        try:
            is_filename = not "\n" in source
        except:
            is_filename = True
        if is_filename:
            with open(source, "rt") as f:
               for l in f:
                   yield l.strip()
        else:
            if source[-1] == "\n":
                 source = source[:-1]
            for l in source.split("\n"):
                 yield l.strip()     



m_entry = re.compile(r"\s*([A-Z_]+)\s+([A-Za-z0-9_]+)\s*;\s*(.+)$")



def make_pxd_entry(line):
    s = ""
    m = m_entry.match(line)
    if m is None:
        return ""
    kwd, par, prototype = m.groups()
    if 'p' in par:
        if kwd == 'EXPORT' and 'x' in par:
            s += "    # PYX <wrap pxi>\n"
        if kwd in ['EXPORT', 'EXPORT_TABLE']:
            s += "    %s\n" % prototype
    return s  
    





def generate_pxd(pxd_out, exp_in, h_name = None, pxd_in = None, nogil = False):
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

    :param exp_in:

        Name of the export file from which external function are copied
        into the .pxd file as a ``cdef extern from`` statement.
        This parameter may be an instance of class ``TableGenerator``.
        The the lines contained in the attribute ``export_file`` 
        of that instance are taken as the lines of the export file    

    :param h_name:
         
         Name of .h file to to be placed into the statement
         ``cdef extern from <h_file_name>``.

    :param pxd_in:
 
         An optional name of a input .pxd file to be copied in the
         output .pxd file in front of the ``cdef`` statement. If
         that parmeter is a string containing a newline character
         ``\n`` then that string is copied directly into the
         output .pxd file.

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
    print('cdef extern from "{0}"{1}:'.format(h_name, s_gil), file = f)
    for line in iter_export_in(exp_in):
        print(make_pxd_entry(line), file = f, end = "")
    if must_open_pxd_out: 
        f.close()


if __name__ == "__main__":
    h_in = os.path.join('..', 'dev', 'c_files', 'clifford12.h')
    generate_pxd(sys.stdout, h_in, nogil=True)





