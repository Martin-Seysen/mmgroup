"""Auxiliary functions for writing tables into C program"""

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



import sys
import re
import os
import warnings
import ast


def pxi_defs(pxifile_name):
    """Return dictionary containing definitions in a .pxi file.

    .pxi files used in the Cython environment contain lines of shape:

    DEF <name> = <literal>

    The function returns a dictionary with entries <name> : <literal>
    """
    f = open(pxifile_name, "rt")
    m_def = re.compile(r"\s*DEF\s+(\w+)\s*=\s*(.+)")
    d = {}
    for l in f:
        try:
            name, value = m_def.match(l).groups()
            d[name] = ast.literal_eval(value.strip())
        except:
           pass
    f.close()
    return d





m_pxd_line = re.compile(r"\s+(\w+)\s+(\w+)\(([^()]*)\)")
m_pxd_arglist = re.compile(r"\s*(\w+)(\s*(\*)?\s*)(\w+)")

def _parse_pxd_line(l):
    m = m_pxd_line.match(l)
    if m:
         type, function, args = m.groups()
         arglist = args.split(",")
         args = []
         for arg in arglist:
             if not arg or arg.isspace():
                 continue
             m_arg = m_pxd_arglist.match(arg)
             if not m_arg:
                 return None
             g = m_arg.groups()
             args.append( (g[0],) + g[2:])
         return (type, function, args) 


m_pxd_enable = re.compile(r"\s*#\s+PYX(.+)")

def _parse_pxd_enable(l):
    m = m_pxd_enable.match(l) 
    if m:
        return True
    else:
        return False



def pxd_to_pxi(pxd_file, module = None, translate = None, nogil = False): 
    """Extract Cython wrappers from prototypes in a .pxd file

    A .pxd file contains prototypes of external C functions and it
    may be included into a .pyx file in order to create a Cython 
    module that uses these C functions.

    This function returns a string containing python wrappers of the
    C functions in the pxd file. This string may be used as a part of 
    an automatically generated .pxi file that can be included directly
    into a .pyx file. Here the C functions in the .pxd file to be 
    included must be sufficiently simple as indicated below.

    The returned string starts with two lines::
 
        cimport cython 
        cimport <module>
    
    Here <module> is given by the parameter ``module``, By default,
    ``module`` is constructed from the name ``pxd_file`` of the .
    pxd file.

    The .pxd file is parsed for lines of shape::

       <return_type> <function>(<type1> <arg1>,  <type2> <arg2> , ...)

    Every such line is converted to a string that codes function,
    which is a Cython wrapper of that function of shape::

       @cython.boundscheck(False)  # Deactivate bounds checking
       @cython.wraparound(False)   # Deactivate negative indexing.
       def <translated_function>(<arg1>, <arg2> , ...):
           cdef <type1> <arg1> = <arg1>_v_
           cdef <type2> <arg2> = <arg2>_v
           cdef <return_type>  _ret
           ret_ = <module>.<function>(<arg1>_v_,  <arg2>_v_ , ...)
           return ret_

    ``<translated_function>`` is the string computed as 
    ``translate(<function>)``, if the argument ``translate`` is given.
    Otherwise is ``<translated_function>`` is equal to ``<function>``.
    
    ``<return_type>`` and ``<type1>, <type2>, ...`` must be valid types 
    known to Cython and C. The types  ``<type1>, <type2>`` used for 
    arguments may  also be pointers. Then pointers are converted to
    memory views e.g::

        uint32_t <function>(uint8_t *a)

    is converted to::

       @cython.boundscheck(False)  # Deactivate bounds checking
       @cython.wraparound(False)   # Deactivate negative indexing.
       def <translated_function1>(a):
           cdef uint8_t[::1] a_v_ = a
           cdef uint32_t ret_
           ret_ = <function>(&a_v_[0])
           return ret_

    The ``<return_type>`` may be void, but not be a pointer. Other 
    types  are not allowed.

    We convert a function only
    if a line containing the string "``# PYX``" (possibly with
    leading blanks) precedes the declaration of a function.
    In the code generator you may 
    use the ``EXPORT`` directive with options ``px`` to enter such
    a line into the source file.

    If ``nogil`` is True, a C function is called with as follows::

       with nogil:
           ret_ = <function>(&a_v_[0])

    Then the C function  ``<function>`` must be e.g. declared as 
    follows::
 
        cdef extern from "<file>.h" nogil:
            int <function> (int a, int *b)
      
    This feature releases the GIL (Global Interpreter Lock) when
    calling the function.
    """
    if not module:
        module = os.path.splitext(os.path.basename(pxd_file))[0]
    s = "cimport cython\n"
    s += "cimport %s\n\n" % module
    if nogil:
        s += "cimport cython\n"
    enable = False
    with open(pxd_file, "rt") as input_file:
        for l in input_file:
            data = _parse_pxd_line(l)
            if not data or not enable:
                enable = _parse_pxd_enable(l)
                continue
            return_type, function, args = data
            t_function = translate(function) if translate else function
            s += "\n"
            s += "@cython.wraparound(False)\n"
            s += "@cython.boundscheck(False)\n"
            s += "def {f}({args}):\n".format( 
              f = t_function, args = ", ".join([a[2] for a in args]))
            c_args = []
 
            has_returnvalue = return_type != "void"
            for argtype, is_ptr, name in args:
                memview = "[::1]"  if is_ptr else ""
                s += "    cdef {argtype}{memview} {name}_v_ = {name}\n".format(
                    argtype=argtype, memview=memview, name=name)
                if is_ptr:
                    c_args.append("&{name}_v_[0]".format(name=name))
                else:
                    c_args.append("{name}_v_".format(name=name))
            assign = ""
            if has_returnvalue:
                s += "    cdef {type} ret_\n".format(type=return_type)
                assign = "ret_ = "
            if nogil:
                s += "    with nogil:\n    "
            s += "    {assign}{m}.{f}({args})\n".format(
                assign = assign, m = module, f = function, 
                args = ", ".join(c_args)
                ) 
            if has_returnvalue:
                s += "    return ret_\n"
            enable = False
    return s




