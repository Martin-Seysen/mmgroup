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

If hundreds of threads use the default seed, it may be useful
to call function ``rand_set_seed(None)`` from time to time.

A seed created by function ``rand_make_seed`` may be initalized
from a volatile source or from a fixed source (which is a 
64-bit integer). Such a seed may be used by one thread only.


Python functions in module ``mmgroup.generators``
.................................................


The following functions should be imported from
module ``mmgroup.generators``.


.. py:function:: .rand_get_seed(seed = None)
    :noindex:
 
    Return seed associated with parameter ``seed`` 

    This function is used for obtaining a suitable python
    object for calling a C function that deals with the 
    random generator. In case ``seed = None`` (default) it 
    returns the default seed object.

.. py:function:: .rand_make_seed(value = None)
    :noindex:

    Create a seed object for the random generator

    The function creates a seed object and returns that object.
    It is intialized with parameter ``value``, which must be
    ``None`` (default) or an unsigned 64-bit integer.

    By default, the seed is initilized from a volatile source. 

    The returned seed object is not thread safe.
 

.. py:function:: .rand_set_seed(seed = None, value = None)
    :noindex:

    Reseed an existing seed object

    The function reseeds an existing seed object ``seed``
    with a ``value``.

    Parameter ``seed`` may be ``None`` (default, referring to 
    the default seed) or an existing seed object created
    by function ``rand_make_seed``.

    Parameter ``value`` is a value for reseeding a seed object
    as in function ``rand_make_seed``. 

    In case ``seed = None`` (default) parameter ``value``
    must also be ``None`` (default).

    If many threads use the default seed then one default
    seed object is created per thread. There is no easy way to 
    determine when a thread has terminated. Here calling
    ``rand_set_seed(None)`` gives all these seed objects
    to the garbage collector, and new default seed objects 
    will be created on demand.


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
