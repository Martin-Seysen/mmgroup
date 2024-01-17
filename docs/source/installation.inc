
Installation on 64-bit Windows
------------------------------

For installing the *mmgroup* package on a 64-bit Windows system,
python 3 must be installed. Then type:

.. code-block:: text

   pip install mmgroup
   pip install pytest
   python -m pytest --pyargs mmgroup -m "not slow"


The last command tests the installation.

Installation on Linux and macOS with a 64-bit x86 CPU
------------------------------------------------------

For installing the *mmgroup* package on a Linux or macOS system,
python 3 must be installed. Then type:


.. code-block:: text

   pip3 install mmgroup
   pip3 install pytest
   python3 -m pytest --pyargs mmgroup -m "not slow"


The last command tests the installation.


Other platforms
---------------

For other platforms the package must be compiled from a
source distribution. Therefore the *cibuildwheel* tool
may be used, see :ref:`cibuildwheel_label`. 

A general description of the build process is given in the
*Guide for developers* of the project documentation.

