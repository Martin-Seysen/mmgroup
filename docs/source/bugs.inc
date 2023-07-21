

In this section we decribe some typical bugs that may occur during the
develpment process and that may be extremely hard to find.



1. Cleaning up a corrupted local git repository
-----------------------------------------------

After a long development session the local git repository might be corrupted
for what reason ever. Here the best cure is to delete everything in the 
repository (exept for subdirectory .git) and to checkout everythinh with git.
If you do not want to do so, you may cleanaup intermediate files as follows.

First open a command shell, switch to the root directory of your local
repository, and make sure that your C compiler works in that shell. Then 
enter the following statements:

 .. code-block:: none

     python cleanup.py -a
     git checkout .
     python setup.py  build_ext --inplace     
     python -m pytest ./src/mmgroup/ -v -s -m "not slow"

Then for creating the documentation you should enter something
like the following (details depending on the operating system):

 .. code-block:: none

     cd docs
     sphinx-build -M html source build  -E -a
     sphinx-build -M latexpdf source build  -E -a
     cd ..

2. Documentation is created correctly on the local host but not on readthedocs
------------------------------------------------------------------------------

Whenever we update the *master* branch in *git* on the *gitbub* server, this 
automatically starts a process on the *readthedocs* server that updates 
the documentation.
In that case the developer should check the documentation on the
*readthedocs* server at

https://mmgroup.readthedocs.io/en/latest/index.html

If this documentation is incorrect then the developer should check the 
logfiles of the last generation of the documentation on *readthedocs*. 
For generating documentation, the *readthedocs* server starts a sequence
of processes, and for each process the output is logged.

There is one type of bug that is difficult to find. If this bug occurs,
then in the output of one of the processes that invokes
``python -m sphinx ...`` there may be an error message  similar to the 
following message:

WARNING: autodoc: failed to import module 'mm_group' from module 'mmgroup'; 
the following exception was raised:

No module named 'mmgroup.mat24'

This message has to do with the mockup of python extensions when running
the *Sphinx* documentation tool on the *readthedocs* server.

We do not use any C compiler on the readthedocs server, so that python
extensions written in *Cython* are not available. When running *Sphinx*
then all python objects to be documented are imported, so that their 
docstrings can be processed. Therefore all Cython extensions must be mocked 
up in *Sphinx*. This means that we tell *Sphinx* that these extensions are 
not available. For details see

https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html#confval-autodoc_mock_imports

The Cython extension ``mat24`` is not mocked up in *Sphinx*, but it
is replaced by class ``mmgroup.dev.mat24.mat24_ref.Mat24`` instead.

So instad of simply importing the extension ``mat24`` we have to do
the following:

 .. code-block:: python

  >>> try:
  >>>    from mmgroup import mat24
  >>> except (ImportError, ModuleNotFoundError):
  >>>    from mmgroup.dev.mat24.mat24_ref import Mat24  
  >>>    mat24 = Mat24  
  

More sophisticated examples of importing objects from module ``mat24``
are given e.g. in modules ``mmgroup.mm_group`` and
``mmgroup.structures.autpl``. 
  
One way to debug such a problem on the local host (and not on the
*readthedocs* server) is as follows. Use ``pip uninstall mmgroup``
(or similar) to uninstall the ``mmgroup package`` Then execute the 
following commands in your shell: 

 .. code-block:: none

     python cleanup.py -a
     git checkout .
     set readthedocs=True
     python setup.py  build_ext --inplace     
     cd docs
     sphinx-build -M html source build  -E -a
     sphinx-build -M latexpdf source build  -E -a
     cd ..
     # Here you should inspect the output of the sphinx-build
     # commands for debugging
     unset readthedocs  # in Windows use 'set readthedocs=' instead  
     # Reverse the effect of  'set readthedocs=True' 
     
If you are lucky then the errors occuring on the *readthedocs* 
server will now also occur on your local host.

.. warning:: 

    After typing that sequence into your shell, your environment
    is probably corrupted; and you should proceed as in subsection
    **Cleaning up a corrupted local git repository**.