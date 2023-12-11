r"""

The python extension ``mmgroup.generators`` provides a random
generator. Its main purpose is the very fast generation of
a vector of integers modulo :math:`p` for the representation
:math:`\rho_p` of the monster group for small numbers :math:`p`.

The internal operation of the random generator is described in
file ``gen_random.c``. Here we describe the python interface
of the random generator. 

This interface is given by a set of functions that can be
imported from python module ``mmgroup.generators``.

A python function dealing with the randon generator requires
a seed. Here a seed is either an object returned by function
``rand_make_seed`` or ``None`` (default).

The default seed is (hopefully) thread save, and it is
initialized from volatile sources such as the time, the process
and the thread id, etc.

A seed created by function ``rand_make_seed`` is initalized
from a fixed source (which is a 64-bit integer). 

Each seed may be used by one thread only. In python a seed is
crreated for each thread whenever needed.


Python functions in module ``mmgroup.generators``
.................................................


The following functions should be imported from
module ``mmgroup.generators``.


.. py:function:: .rand_get_seed()
    :noindex:
 
    Return the seed associated with the current thread 

    This function is used for obtaining a suitable python
    object for calling a C function that deals with the 
    random generator. The user must not modify that seed!


.. py:function:: .rand_make_seed(valu )
    :noindex:

    Create a deterministic seed object for the random generator

    The function creates a seed object and returns that object.
    It is intialized with parameter ``value``, which must be
    an unsigned 64-bit integer.
 

.. py:function:: .rand_bytes_modp(p, num_bytes, seed = None)
    :noindex:

    Return a random array of integers mod p
  
    The function returns a one-dimensional numpy array 
    of uniform distributed random integers ``x`` with 
    ``0 <= x < p``. That array has length ``num_bytes``
    and has ``dtype = numpy.uint8``. 

    Parameter ``seed`` must be either ``None`` (default, 
    referring to the default seed) or a seed object
    created by function ``rand_make_seed``.


.. py:function:: .rand_fill_bytes_modp(p, array_bytes, seed = None)
    :noindex:

    Fill an existing array with random integers mod p
  
    The function fills the one-dimensional numpy array ``array_bytes``
    with  uniform distributed random integers ``x`` with ``0 <= x < p``. 
    That array must have ``dtype = numpy.uint8``. 

    Parameter ``seed`` must be either ``None`` (default, 
    referring to the default seed) or a seed object
    created by function ``rand_make_seed``.


.. py:function:: .rand_gen_modp(p, seed = None)
    :noindex:

    Return random integer ``x`` with ``0 <= x < p``.

    Here 1 <= p <= 2**32 must hold.

    Parameter ``seed`` must be either ``None`` (default, 
    referring to the default seed) or a seed object
    created by function ``rand_make_seed``.


"""
