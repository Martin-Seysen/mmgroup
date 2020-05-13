

The current version of the ``mmgroup`` package is a source distribution that has 
been tested on a 64-bit Windows platform only. It runs with python 3.6 or higher. 
The source distribution can be downloaded from

https://github.com/Martin-Seysen/mmgroup . 

The distribution contains a number of extensions written
in ``C`` which have to be built before use. At present there is no binary 
distribution for Windows.

Dependencies
------------

Before you can use this package or build its extensions you should install 
the following python packages:

.. table:: External Python packages required
  :widths: 30 70


  ======================== =========================================================
  Package                  Purpose 
  ======================== =========================================================
  ``cython``               Development: integrating ``C`` programs into the 
                           ``mmgroup`` package
  ``numpy``                Runtime: Most arrays used by the ``mmgroup`` package 
                           are ``numpy`` arrays
  ``pytest``               Testing: basic package used for testing
  ``scipy``                Testing: statistical tests
  ``setuptools``           Development: basic package used for setup and building
                           extensions
  ``sphinx``               Documentation: basic package used for documentation
  ``sphinxcontrib-bibtex`` Documentation: bibliography in BibTeX style
  ======================== =========================================================

Packages used for the purpose of documentation are required only if you want to 
rebuild the documentation. If you want to rebuild the documentation you should
also install the following programs:


.. table:: External programs required
  :widths: 20 30 50


  ========= ============== ===========================================================
  Program   Purpose        Location 
  ========= ============== ===========================================================
  miktex    Documentation  https://miktex.org/
  Perl      Documentation  https://www.perl.org/get.html
  ========= ============== ===========================================================



Installing the extension
------------------------

To install the required extension, go to root directory of the distribution.
This is the directory containing the files ``setup.py`` and ``README.rst``.
From there run the following command::

    python setup.py install -f

For testing the installation, run the following command::

    pytest src/mmgroup/ -v -m "not slow"









Installing a C compiler for cython in Windows 
---------------------------------------------

The bad news for Windows developers is that there is no pre-installed 
``C`` compiler on a standard Windows system. However, the ``cygwin`` 
package requires a ``C`` compiler. Here in principle, the user has the 
choice between the following two compilers:

   *  ``MinGW-w64``

   *  ``MSVC``

The user has to install a ``C`` compiler so that it cooperates with 
``cython``. 
That installation process is out of the scope of this document.

For installing ``MinGW``, one might start looking at
https://cython.readthedocs.io/en/latest/src/tutorial/appendix.html.

For installing ``MSVC``, one might start looking at 
https://wiki.python.org/moin/WindowsCompilers

The current ``setup.py`` supports ``MinGW`` for 64-bit Windows only.

The 64-bit Anaconda distribution of Python 3.7 appears to
contain a suitable ``MinGW`` version.

