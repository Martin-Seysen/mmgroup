

================================
The mmgroup guide for developers
================================


Introduction
============

This *guide for developers* explains some selected topics of the
``mmgroup`` project which are relevant for developers who want to make 
changes or bugfixes. It also explains some mathematical aspects of the 
implementation which are not covered by :cite:`Seysen20`.

The *guide* is far from complete and we hope that we can provide 
more information in future versions of it.

Directory structure
===================

.. include:: guide_directory.inc


.. _install_from_source_label:

Installation from a source distribution
=======================================


.. include:: installation_from_source.inc


The build process
=================

.. include:: guide_build_process.inc


Some mathematical aspects of the implementation
===============================================

.. include:: guide_mathematical.inc


.. include:: guide_mathematical_2.inc


.. _code-generation-label:

Code generation
===============

Warning!

At present the process of generating C code is under construction.

We plan to switch to the ``meson`` build system. Therefore a much more
declarative style is required for describing an buld operations.
Thus the description here is pertty much outdated!


In this section we describe the most important functions and classes
used for the automatic generation of C code.
  

.. automodule:: mmgroup.generate_c.make_c_tables_doc


Classes and functions provided by the code generator
----------------------------------------------------
 

.. autoclass:: mmgroup.generate_c.TableGenerator
   :members: generate


.. autoclass:: mmgroup.generate_c.UserDirective

.. autoclass:: mmgroup.generate_c.UserFormat

          
.. autofunction::   mmgroup.generate_c.c_snippet  

.. autofunction::   mmgroup.generate_c.format_item

.. autofunction::   mmgroup.generate_c.make_doc

.. autofunction::   mmgroup.generate_c.make_table

.. autofunction::   mmgroup.generate_c.prepend_blanks
   
.. autofunction::   mmgroup.generate_c.pxd_to_pxi
  
.. autofunction::   mmgroup.generate_c.generate_pxd

    

How the code generator is used
------------------------------

.. include:: guide_use_code_generation.inc



Description of some typical bugs
================================

.. include:: bugs.inc

