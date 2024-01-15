

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



Some mathematical aspects of the implementation
===============================================

.. include:: guide_mathematical.inc


.. include:: guide_mathematical_2.inc



.. _install_from_source_label:

Installation from a source distribution
=======================================


.. include:: installation_from_source.inc



.. _build-process-label:


The build process
=================


Warning! At present the build process is under construction.


We plan to switch to the ``Meson`` build system. The reason for this
is that invoking ``setup.py`` from the command line is deprecated.
``Meson`` requires a much more declarative style for describing the 
build operations.
Thus the description here is pretty much outdated!


.. include:: guide_build_process.inc



.. _code-generation-label:

Code generation
===============

Warning! At present the code generation process is under construction.

We plan to switch to the ``Meson`` build system. Therefore a much more
declarative style is required for describing the code generation.
Also, a strict separation between input and output directories is
required.
Thus the description here is pretty much outdated!



.. include:: guide_code_generation.inc



How the code generator is used
------------------------------

.. include:: guide_use_code_generation.inc



Description of some typical bugs
================================

.. include:: bugs.inc

