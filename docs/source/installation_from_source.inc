

The current version of the ``mmgroup`` package is a source distribution 
that has been tested on a bit Windows, Linux and macOS with a 
64-bit x85 CPU. It runs with python 3.8 or higher. The sources of the
project can be downloaded from

https://github.com/Martin-Seysen/mmgroup . 


We use the python package ``cibuildwheel`` to build several wheels on
these operating systems. The GitHub repository of the project contains
actions to build the corrsponding python wheels. The program code for 
the GitHub actions is stored in subdirectory ``.github/workflows``.
More details about the build process are given in
section :ref:`build-process-label`.

The *mmmgroup* package contains a number of extensions written
in ``C`` which have to be built before use. This will be discussed
in section :ref:`code-generation-label`.


.. _cibuildwheel_label:

Building the ``mmgroup`` package with cibuildwheel
---------------------------------------------------

The easiest way to build the ``mmgroup`` package is to clone the sources
from the github repository with *git*, and to build a python wheel
with *cibuildwheel*. Therefore, the cibuildwheel tool must support
your target platform. In the sequel we assume that python and the git
revision control system are installed on your platform.
We describe the build process on a unix-like platform.

Depending on your platform, you might have to install additional tools.
For details we refer to the cibuildwheel docmentation
https://cibuildwheel.readthedocs.io/en/stable/

Any previously installed version of the ``mmgroup`` package must
be uninstalled before building a new version of that package!

E.g. for building the ``mmgroup`` package for python 3.12 on Linux 
with a 64-bit x86 CPU, open a shell and type:

.. code-block:: bash

    pip3 install cibuildwheel
    git clone https://github.com/Martin-Seysen/mmgroup
    cd mmgroup
    git pull origin master
    git checkout .
    python3 -m cibuildwheel --output-dir wheelhouse --only cp312-manylinux_x86_64

The last command builds the python wheel. For building the wheel on
a different platform, the argument 'cp312-manylinux_x86_64' in the
last command should be replaced by a *build identifier* that
describes the requested python version and target platform.
A list of valid build identifiers is given in Section 
*'Options/build selection'* of the cibuildwheel documentation, see

https://cibuildwheel.readthedocs.io/en/stable/options/#build-skip

After building the python wheel, switch to subdirectory *wheelhouse*,
and check the file name of the wheel that has been built:

.. code-block:: bash

    cd wheelhouse
    ls

Then you may install the wheel with pip, e.g.:  

.. code-block:: bash

    pip3 install mmgroup-0.0.13-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl


Here the file name after the *pip3 install* command must be changed 
to the file name of the wheel just generated. 

You should test your installation of the ``mmgroup`` package before
using it:

.. code-block:: bash

   pip3 install pytest
   python3 -m pytest --pyargs mmgroup -m "not slow"


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
  ``auditwheel``           Development (Linux and macOS only): patching shared
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


Remarks
.......

If you have installed any version of the ``mmgroup`` package (e.g. with
the ``pip`` tool) then you **must** uninstall that package before you can 
build a new version of the ``mmgroup`` package from the source files.



Installing a C compiler for cython in Windows 
---------------------------------------------

The bad news for Windows developers is that there is no pre-installed 
``C`` compiler on a standard Windows system. However, the ``cython`` 
package requires the ``C`` compiler ``MSVC``.   


The user has to install a ``C`` compiler so that it cooperates with 
``cython``. 
That installation process is out of the scope of this document.

For installing ``MSVC``, one might start looking at 
https://wiki.python.org/moin/WindowsCompilers

The author has installed the ``MSVC`` compiler with the Microsoft
*Build Tools for Visual Studio* from:

https://visualstudio.microsoft.com/thank-you-downloading-visual-studio/?sku=BuildTools&rel=16 ,

following the instructions in

https://www.scivision.dev/python-windows-visual-c-14-required/  .

Before typing ``python setup.py bdist_wheel`` in a Windows command 
line the author had to type::

  "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Auxiliary\Build\vcvars64.bat"

Here the path my be different on the user's Windows system.



    

