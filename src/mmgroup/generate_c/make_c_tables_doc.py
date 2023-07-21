r"""The C code generator generates code from a source file automatically.


    Class ``TableGenerator`` can be used to generate  a 
    .c and a .h file from *source* file. Here such a source file usually
    has the extension .ske.
 
    Basic idea: We copy a *source* file to the .c file. The source 
    contains directives to enter precomputed tables, code  snippets 
    or constants into the .c file, and prototypes into the .h file.

    We copy the content of the *source* file directly to the .c
    file.   In the *source* file we also parse certain directives 
    and we replace them by automatically generated code snippets 
    in the .c file and in the .h file. 

    The arguments given to the constructor of class ``TableGenerator``
    specify the meaning of the directives.

    Overview
    --------   

    In the *source* file we parse directives of the shape::

      // %%<keyword>  <arg1>, <arg2>, ...

    Each directive executes certain function associated to the
    ``<keyword>``. There are built-in and user-defined directives.
    Built-in directives are listed in section *Built in directives*.
    User-defined directives are  passed to the constructor in a
    dictionary ``directives`` with entries::

            <keyword> : <directive_function>. 
    
    If a directive is found in the source file we call::

       <directive_function>(*<all evaluated arguments>), 

    and the string returned by that function is merged into the
    generated  .c  file. The arguments of a directive must be given 
    in python syntax and they  are evaluated to python objects using 
    some safe version of the python ``eval()`` function. Basically, 
    python identifiers in an argument are evaluated as follows:

    In the constructor, parameter ``tables`` is a dictionary of local
    variables, which will be passed to the python ``eval()`` function 
    for evaluating an argument. See section *Arguments of directives* for
    a more detailed description of this process.
      
    You may also consider class ``UserDirective`` 
    for constructing specific directive functions.


    Arguments of directives
    -----------------------
     
    The comma-separated list of arguments of a built-in or user-defined 
    directive or function must be a valid python expression. This 
    expression is evaluated by a safe version of the python ``eval()`` 
    function,  with the local variables given by the dictionary 
    ``names``. 
    Here ``names`` is  an updated version of the dictionary ``tables`` 
    of local variables  passed in the constructor of this class as
    described below. This means  that each identifier found in the list  
    of arguments is looked up in that dictionary, and that the identifier 
    is evaluated  to the corresponding value in that dictionary.

    The dynamic dictionary ``names`` is  initialized with dictionary  
    ``tables``. It may be updated  with the names of the tables generated 
    by the built-in directive ``TABLE``, as indicated in the
    description of the  built-in directive ``USE_TABLE``.

    There are a few more predefined entries in dictionary ``names``:

    ==========  ==============================================================
    name        value
    ==========  ==============================================================
    ``TABLES``  The original dictionary ``tables`` passed  in the  constructor
    ``NAMES``   updated version of dictionary ``tables`` as described above
    ==========  ==============================================================

    Evaluating a python expression with the standard function ``eval()`` 
    is known to be unsafe. For evaluating arguments we use a modified 
    version of the ``eval()`` function. That modified function does not
    allow assignment, and identifiers starting with an underscore ``_`` 
    are illegal. See function ``EvalNodeVisitor`` in module
    ``mmgroup.generate_c.generate_functions``. for details. 
    
    There is a limited number of built-in python functions (specified 
    in directory ``safe_locals`` in the same module) that are
    legal in python expressions for arguments::

       abs, int, len, max, min, str, pow, range, zip 

    An example with user-defined  directives is given in class 
    ``testcode.TestGen``. We recommend  to follow the conventions in that 
    class for passing user-defined directives to an instance of class 
    ``TableGenerator``.   


    Built-in directives
    -------------------

     * ``COMMENT``,  deprecated!
        Indicates that a comment follows. This should be used only
        for comments relevant for the user, not for internal details.
        No action in this class, but other functions may use this
        directive for generating user documentation.  

     * ``ELSE [IF]?`` <expression>
        *else* clause for IF directive, see ``IF``

     * ``END`` <directive>
        Closes a block that has been started by one of the directives
        ``FOR`` or ``IF``. A ``FOR`` block must be closed with ``END FOR``, 
        an ``IF`` block must be closed with ``END IF``, etc.

     * ``ERROR``  <message>
         TODO: yet to be implemented and documented!!!!


     * ``EXPORT``   (arguments see below)
        Then the following line with an appended semicolon is copied 
        into the generated .h file. This is useful for 
        creating prototypes.
        The ``EXPORT`` directive may precede a line with a ``USE_TABLE`` 
        directive.
        Then the table name captured by the next ``USE_TABLE``  statement 
        is  exported to the generated .h file. 
        The ``EXPORT`` directive may precede a line with a  ``DEFINE``  function.
        Then the #define statement generated by the next DEFINE function
        is written to the generated .h file.

        An optional argument ``p`` means that the exported line will also
        be included in a .pxd file, see method ``generate_pxd()``.

        An optional argument ``px`` or ``xp`` means that the same as 
        argument ``p``. Then in addition we also write a comment 
        ``# PYX <wrap pxi>`` into the .pxd file generated by method 
        ``generate_pxd()``. If function ``mmgroup.generate_c.pxd_to_pyx``
        reads such a comment in a .pxd file, it will write a simple wrapper
        for the exported function into a .pyx file.


     * ``EXPORT_KWD`` <keyword>
        The directive:: 

          // %% EXPORT_KWD  FOO_API

        places the keyword FOO_API before any function or variable 
        which is exported with the ``EXPORT`` or ``EXPORT_TABLE``
        directive.
        This is useful on a MS Windows system for creating a dll.
        Then the header should contain a definition similar to::

          #define FOO_API __declspec(dllexport)

        This declares these functions or variables to be exported by the
        dll. For background, see e.g. 
        https://gcc.gnu.org/wiki/Visibility.

     * ``EXPORT_TABLE`` (no arguments)
        Parses the next input line for the C-name of a table, and 
        exports this name into the generate .h header file.
        This is equivalent to two subsequent directives::

           EXPORT
           USE_TABLE

        There are, however, subtle differences, see directive 
        ``EXPORT_KWD``.

     * ``FOR``  <variable> ``in`` <parameter_list>
        Generate a sequence of blocks of statements, one block for 
        each parameter in the parameter_list. A sequence of arbitrary
        lines follows, terminated by an ``END FOR`` directive. The
        parameter list is evaluated to a python object which
        should be iterable.

        Nesting of FOR and similar directives is possible. 
  
        Examples::   

         // %%FOR x in range(1,3)
         i%{x} += 2*%{x};
         // %%END FOR

        evaluates to::

         i1 += 2*1;
         i2 += 2*2;
            

     * ``GEN`` <ch>
        Directs output to .c or to .h file or to both, depending on
        whether the letter c or h is present in the first argument.

     * ``IF`` <expression>
        Conditional code generation depending on the boolean value
        of <expression>.  Syntax is::

          // %%IF <expression>
          <statements>
          // %%ELSE IF <expression>
          <statements>
          // %%ELSE 
          <statements>
          // %%END IF

        with one or more optional ``ELSE IF`` clauses and at most one
        final ``ELSE`` clause. Nested ``IF`` and ``FOR`` blocks are
        possible.


     * ``INCLUDE_HEADERS``

        This directive is legal in a header file with extension ``.h``
        only. A typical header file may look like this::

          #ifndef  _THIS_HAEDER_HAS_ALREADY_BEEN_PROCESSED_
          #define  _THIS_HAEDER_HAS_ALREADY_BEEN_PROCESSED_

          // Some prototypes for C functions will be inserted here

          #endif

       The code generator may insert automatically generated prototypes
       for C functions into that header. You can tell the code generator
       where to insert these prototypes as follows::

          #ifndef  _THIS_HAEDER_HAS_ALREADY_BEEN_PROCESSED_
          #define  _THIS_HAEDER_HAS_ALREADY_BEEN_PROCESSED_

          // %%INCLUDE_HEADERS

          #endif
        

     * ``JOIN`` <infix>, <suffix>    
        A ``FOR`` directive must follow a ``JOIN`` directive. The 
        expressions generated by the following ``FOR`` directive
        are joined with the <infix> (which is usually an
        operator), and the <suffix> is appended to the
        whole generated expression. E.g.::

          a += 
          // %%JOIN*  " +", ";"
          // %%FOR* i in range(1,4)
            %{int:2*i} * table[%{i}]
          // %%END FOR

        evaluates to::
 
          a += 
            2 * table[1] + 
            4 * table[2] + 
            6 * table[3];

     * ``PY_DOCSTR`` <string>, <format>
         This directive is no longer supported       

     * ``PYX`` <string>
         This directive is no longer supported       


        into the .pxd file, when generating a .pxd file. This has
        no effect on the generated .c or .h file. It may be used for
        automatically generating a .pyx file from a .pxd file.

     * ``TABLE``  <table>, <format>
        Enters all entries of the python array <table> into the C code 
        as an array of constants with an optional <format>. Only the 
        data of the array are written. Feasible <format> strings are 
        given in function format_number().

     * ``USE_TABLE``   (no arguments)
        Parses the next input line for the C-name of a table and makes 
        that C name available for user-defined directives. 

        If a table is coded  with  the ``TABLE`` directive and
        preceded by a ``USE_TABLE`` or ``EXPORT_TABLE`` directive then
        the  name of that table will be added to the dictionary ``names``.
        The key with the python name of that table will have the
        C name of that table as its value.

        Consider the following example::

         1:    // %%USE_TABLE
         2:    int  C_table[]  = {
         3:       // %%TABLE python_table, uint32
         4:    } 

        Here ``'python_table'`` must be a key in the dictionary ``tables``
        passed to the constructor of class ``TableGenerator``. The
        value for that key should be a list of 32-bit integers.
        Then the directive in line 3 creates the list of integer
        constants given by the list ``tables[python_table]``.

        The ``USE_TABLE`` directive in line 1 has the following effect.
        In the sequel the entry ``names['python_table']`` has the value
        ``'C_table'``, which is the C name of the table that has been
        generated. So we may use that entry in subsequent directives
        or via string formatting with curly braces.  
      
        The entries of table ``python_table`` are still available in  
        the form ``python_table[i]`` as before. ``python_table[:]`` 
        evaluates to  the whole table as before.

     * ``WITH``  <variable> = <value>
        Temporarily assign  <value> to <variable>. A sequence of 
        arbitrary lines follows, terminated by an ``END WITH`` 
        directive. Both, <variable> and <value>, may be tuples of 
        equal length; a nested tuple <variable> is illegal. 
        The <variable>, or the variables contained in the tuple
        <variable>, are temporarily added to the dictionary ``names`` 
        as keys, to that they can be used via string formatting.

        Nesting of WITH and other directives is possible. 
  
        Example::   

          // %WITH square, cube = int((x+1)**2), int((x+1)**3)
          i += %{cube} + %{square} + %{x};
          // %%END WITH

        Assuming ``x`` has value ``9``, this evaluates to::

          i += 1000 + 100 + 9;
         
            



    String formatting
    -----------------

    A line ``l`` in the source file may contain a string ``s``  with 
    balanced curly  braces and preceded by a '%' character.
    Then the standard ``.format()`` method is applied  to  the string 
    ``s[1:]`` obtained from  ``s`` by dropping the initial '%' 
    character. The  dictionary ``names`` is given to the formatting
    method as given as the list of keyword arguments. 

    E.g. in case ``names['x'] = y``, the expressions ``%{x.z}`` and 
    ``%{x[z]}``  refer  to ``y.z`` and ``y[z]``. Of course, the objects 
    ``y.z`` or ``y[z]`` must exist and evaluate to strings which are 
    meaningful in the generated C file.

    For formatting, the dictionary ``names`` is an updated version of 
    the dictionary ``tables`` passed in the constructor, as described 
    in section *Arguments of directives*.
                 
    We consider this to be a reasonable compromise between the standard 
    C syntax and the very powerful pythonic ``.format`` method for 
    strings. 


    Using functions via string formatting
    -------------------------------------

    User-defined functions can be invoked inside C code via the the
    string formatting operator ``%{function_name:arg1,arg2,...}``.
    Then the function with name ``function_name`` is called with the
    given arguments, and the result of the function is substituted 
    in the C code. 

    Such a user-defined function may be created with class 
    ``UserFormat``. The simplest way  to create such a  user-defined 
    function from a function ``f`` is to create an entry::

        "function_name" : UserFormat(f)

    in the dictionary ``tables`` passed to the constructor. Then
    ``%{function_name:arg1,arg2,...}`` evaluates to 
    ``str(f(arg1,arg2,...))``.
    Here arguments ``arg1, arg2`` must be given in python syntax and 
    they are processed in the same way as the arguments of a directive 
    described in section *Arguments of directives*. Then you may also 
    write ``function_name(arg1,..)`` inside an argument of a directive 
    or  a user-defined function; this expression evaluates to 
    ``f(arg1,...)``.

    Class ``UserFormat`` also contains options to control the 
    conversion of the result of function f to a string.

    There are a few built-in string formatting functions, as given
    by dictionary ``built_in_formats`` in module 
    ``generate_functions``. 

    There are a few predefined functions for string formatting:

      * ``%{int:x}``    returns ``x`` as a decimal integer 
      * ``%{hex:x}``    returns ``x`` as a hexadecimal integer
      * ``%{len:l}``    returns length of list ``l`` as a decimal int     
      * ``%{str:x}``    returns ``str(x)``
      * ``%{join:s,l}`` returns ``s.join(l)`` for string ``s`` and list ``l``

    See dictionary ``built_in_formats`` in module 
    ``mmgroup.generate_c.generate_functions`` for details.

    So, assuming that keys ``'a'`` and ``'b'`` have values ``1`` and ``10`` 
    in the dictionary ``tables`` passed to the code generator, we may 
    write e.g.::

       y%{a} += %{int:3*b+1};

    in the source file. This evaluates to::

       y1 += 31;
        

    Quiet form of a directive
    -------------------------

    For directives such as ``IF``, ``FOR``, ... there are also quiet 
    variants ``IF*``, ``FOR*``, ..., which perform the same operations, 
    but without writing any any comments about the directive into the 
    output file. E.g. in::
         
        // %%FOR i in range(50)
        printf("%d\n", i);
        // %%IF* i % 7 == 0 or i % 10 == 7
        printf("seven\n");
        // %%END IF
        // %%END FOR

    you probably do not want any messy comments about the IF directive
    in your generated code.

    Historical remark
    -----------------

    In versions up to 0.4 string formatting operator was ``{xxx}``
    instead of ``%{xxx}``. This had the unpleasent effect that
    valid C expressions have a different meaning in the input and
    and the output of the code generation process. This situation
    became unbearable after intoducing doxygen for documentation
    of C code. Then expressions like ``\f$ x_{1,1} \f$`` could not
    be coded in a reasonable way the old version.

"""

