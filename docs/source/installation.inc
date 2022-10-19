
Installation on 64-bit Windows
------------------------------

For installing the *mmgroup* package on a 64-bit Windows system,
python 3 must be installed. Then type:

.. code-block:: text

   pip install mmgroup
   pip install pytest
   python -m pytest --pyargs -Wignore mmgroup -m "not slow"


The last command tests the installation.
32-bit Windows systems are not supported.

Installation on Linux and macOS
-------------------------------

For installing the *mmgroup* package on a Linux or macOS system,
python 3 must be installed. Then type:


.. code-block:: text

   pip3 install mmgroup
   pip3 install pytest
   python3 -m pytest --pyargs -Wignore mmgroup -m "not slow"


The last command tests the installation.


Other operating systems
-----------------------

For other operating systems, the package must be compiled from a
source distribution. Some hints for doing so are given in the
*guide for developers*.
