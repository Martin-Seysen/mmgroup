

The directory is organized as follows:

  * ``[Root directory]``
      Contains setup.py and configuration files 

    * ``docs``
       * ``sources``
          Source files for documentation with Sphinx 

    * ``src``
       * ``mmgroup``
         The directory containing the package, and 
         also extensions and shared libraries
          
         * ``dev``     
           Contains stuff needed by developers, including scripts and
           source file for generating C code automatically.  
          
         * ``generate_c``           
           The code generator (for generating C code automatically)


         * ``structures``
           Some basic structures used by the modules in ``mmgroup`` 

         * ``tests``
           Test scripts and some algebraic structures used for testing 
           with ``pytest``

             
Directory  ``src/mmgroup/dev`` has the following subdirectories: 

   * ``c_files``
     Contains automatically generated C code

   * ``clifford12``
     Scripts for generating C code related to the subgroup
     :math:`2^{1+24}.\mbox{Co}_1` of the monster and also to
     a Clifford group resembling that subgroup. 
 
   * ``generators``
     Scripts for generating tables related to the monomial operation
     of the element :math:`\xi` of :math:`\mathbb{M}`, and also 
     for generating C code related to certain subgroups of
     :math:`\mathbb{M}`.

   * ``hadamard``
     Scripts for generating C code related to Hadamard matrices

   * ``mat24``
     Scripts for generating C code related to the Mathieu group
     :math:`\mbox{Mat}_{24}`

   * ``mm_basics``
     Scripts for generating C code related to the representation
     :math:`\rho_p`  of :math:`\mathbb{M}`

   * ``mm_op``
     Scripts for generating C code related to the operation of 
     :math:`\mathbb{M}` on :math:`\rho_p`  for specific values
     :math:`p`

   * ``mm_reduce``
     Scripts for generating C code related to the word shortening
     algorithm and the computation of the order in  
     :math:`\mathbb{M}`.

   * ``pxd_files``
     Contains files with extension ``.pyx, .pxd, .pxi`` used by 
     ``Cython``, which are automatically generated or copied from
     a different source location 

 

Directory  ``src/mmgroup/tests`` has the following subdirectories: 

   * ``groups``
     Slow python implementations of certain groups, including a
     preimage of :math:`\mathbb{M}`, for testing.

   * ``spaces``
     Slow python implementations of certain vector spaces, 
     including the space  :math:`\rho_p`, for testing.

   * ``test_xxx``
     Test scripts (to be executed with ``pytest``) for testing
     module ``xxx``.
         