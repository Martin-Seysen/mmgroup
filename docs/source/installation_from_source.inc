

The current version of the ``mmgroup`` package is a source distribution that has 
been tested on a 64-bit Windows, Linux and macOS. It runs with python 3.7 or 
higher. The sources of the project can be downloaded from

https://github.com/Martin-Seysen/mmgroup . 

The distribution contains a number of extensions written
in ``C`` which have to be built before use. 

Dependencies
------------

Before you can use this source distribution or build its extensions you 
should install the following python packages:

.. table:: External Python packages required
  :widths: 30 70


  ======================== =========================================================
  Package                  Purpose 
  ======================== =========================================================
  ``cython``               Development: integrating ``C`` programs into the 
                           ``mmgroup`` package
  ``numpy``                Runtime: Most arrays used by the ``mmgroup`` package 
                           are ``numpy`` arrays
  ``wheel``                Distribution: package for generating a Python wheel
  ``pytest``               Testing: basic package used for testing
  ``regex``                Development: parsing balanced curly braces
  ``setuptools``           Development: basic package used for setup and building
                           extensions
  ``patchelf``             Development (Linux and macOS only): patching shared
                           libraries
  ``sphinx``               Documentation: basic package used for documentation
  ``sphinx-rtd-theme``     Documentation: 'theme' to be used by sphinx 
  ``sphinxcontrib-bibtex`` Documentation: bibliography in BibTeX style
  ``sphinxcontrib-tikz``   Documentation: link between ``doxygen`` and ``sphinx``
  ``breathe``              Documentation: link between ``TikZ`` and ``sphinx``
  ======================== =========================================================

Packages used for the purpose of documentation are required only if you want to 
rebuild the documentation. If you want to rebuild the documentation you should
also install the following programs:


.. table:: External programs required
  :widths: 20 30 50


  ========== ============== ======================================================
  Program    Purpose        Location 
  ========== ============== ======================================================
  miktex     Documentation  https://miktex.org/
  Perl       Documentation  https://www.perl.org/get.html
  doxygen    Documentation  https://www.doxygen.nl/download.html
             of C files
  Ghostcript Documentation: https://ghostscript.com/releases/gsdnld.html
             using TikZ in    
             Sphinx
  ========== ============== ======================================================



Installing the package
----------------------

To build the required package on your local computer, go to the root directory 
of the distribution.
This is the directory containing the files ``setup.py`` and ``README.rst``.
From there run the following commands in a shell:


.. code-block::

   python -m pip install -r requirements.txt
   python setup.py bdist_wheel

In Linux or macOS you'll have to change the ``python`` command to ``python3``.

For testing the installation, run the following command:

.. code-block::

   python -m pytest ./src/mmgroup/ -Wignore -v -s -m "not slow"



Distributing a *wheel* is the standard way to distribute a python
package, see e.g. 

https://packaging.python.org/guides/distributing-packages-using-setuptools/#wheels

The current version supports building wheels for 64-bit Windows, Linux and maccOS.

Remarks
.......

If you have installed any version of the ``mmgroup`` package (e.g. with
the ``pip`` tool) then you **must** uninstall that package before you can 
build a new version of the ``mmgroup`` package from the source files.

We recommend to run

.. code-block::

   python cleanup.py -pcx

in the root directory of the distribution for cleaning up intermediate
files before building the package from the source files.


Installing a C compiler for cython in Windows 
---------------------------------------------

The bad news for Windows developers is that there is no pre-installed 
``C`` compiler on a standard Windows system. However, the ``cython`` 
package requires a ``C`` compiler. Here in principle, the user has the 
choice between the following two compilers:

   *  ``MSVC``   

   *  ``MinGW-w64``

The user has to install a ``C`` compiler so that it cooperates with 
``cython``. 
That installation process is out of the scope of this document.

For installing ``MSVC``, one might start looking at 
https://wiki.python.org/moin/WindowsCompilers

For installing ``MinGW``, one might start looking at
https://cython.readthedocs.io/en/latest/src/tutorial/appendix.html.

The current ``setup.py`` supports ``MSVC`` and ``MinGW`` for 64-bit 
Windows. According to the last URL the ``MinGW`` compiler works with 
Python versions up to 3.4 only, so that it can no longer be recommended. 

The author has installed the ``MSVC`` compiler with the Microsoft
*Build Tools for Visual Studio* from:

https://visualstudio.microsoft.com/thank-you-downloading-visual-studio/?sku=BuildTools&rel=16 ,

following the instructions in

https://www.scivision.dev/python-windows-visual-c-14-required/  .

Before typing ``python setup.py bdist_wheel`` in a Windows command 
line the author had to type::


  "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Auxiliary\Build\vcvars64.bat"

Here the path my be different on the user's Windows system.



    
Remarks
....... 

You may also use the python package ``cibuildwheel`` to build several wheels on
differnt operating systems. The GitHub repository of the project contains
actions to build several wheelson differnt opating systems. 

