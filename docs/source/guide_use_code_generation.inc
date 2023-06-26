Warning!

At present the process of generating C code is under construction.

We plan to switch to the ``meson`` build system. Therefore a much more
declarative style is required for describing an buld operations.
Thus the description here is pertty much outdated!


For generating C code, an instance of class  ``TableGenerator`` in module
``mmgroup.generate_c`` must be generated. That instance takes two 
arguments ``tables`` and ``directives`` as described in section 
:ref:`code-generation-label`.

There are many tables and directives needed for different purposes. By 
convention we provide a variety of classes, where each class provides 
two attributes ``tables`` and ``directives`` providing dictionaries of 
tables and directives, respectively. These dictionaries  are suitable
for creating instances of classes ``TableGenerator``. Such classes
are called *table-providing classes*. Attributes ``tables`` and 
``directives`` may also be implemented as properties.

Python scripts with names ``codegen_xxx.py`` in the root directory of
the source distribution create instances of class ``TableGenerator``
that will generate  C programs and also the corresponding header files. 
These python scripts are run as subprocesses of the ``build_ext`` 
command in the ``setup.py`` script.

Each instance of class ``TableGenerator`` may take the union of the
directories obtained from methods ``tables()`` and ``directives()`` 
of several table-providing classes. A simple example for a 
table-providing class is the class ``Lsbit24Function``
in module ``mmgroup.dev.mat24.mat24_aux``.


