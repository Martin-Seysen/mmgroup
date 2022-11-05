r"""The automatically generated file ``mat24_functions.c`` contains the C
code for the functions exported by the python extension ``mmgroup.mat24``.

These functions are documented in file ``mat24_functions.c``. Here 
just give an overview of the functionality of that module.

File  ``mat24_functions.c`` has been generated from the source file 
``mat24_functions.ske`` using the code generator in module 
``mmgroup.generate_c``. 

The functions in file ``mat24_functions.c`` perform basic computations
in the Golay code, in its cocode, and in the Mathieu group ``Mat24``. 
They also deal with  Parker loop ``Pl`` and with its automorphism group
``AutPl``. A more comfortable python interface to these objects is 
provided by the python classes ``GCode``, ``Cocode``, ``PLoop``, and 
``AutPL`` in module ``mmgroup``.

All C functions in file ``mat24_functions.c`` start with the prefix
``mat24_``. These functions are also called from C functions in other
modules. Therefore we store the binary code for these functions in
a shared library. For some of these functions there are also
macros (defined with ``#define``), starting with ``mat24_def_``.

There is a one-to-one correspondence between the functions in
``mat24_functions.c`` and the function exported from the python
extension ``mmgroup.mat24``, see subsection
*Mapping C functions to python functions* for details.

The python class ``Mat24`` in module ``mmgroup.dev.mat24.mat24_ref``
contains pure python implementations of most functions of the 
``mmgroup.mat24`` extension as class methods. This class is used for 
testing and as a substitute for the ``mmgroup.mat24`` extension in
an early stage of the build process.

In the following subsections the term *C functions* refers to the
C functions in file ``mat24_functions.c``. In the documentation, the 
names of the C functions  are given without the ``mat24_`` prefix.
The term *API reference* means the main document
*The mmgroup API reference* of this project.

The Golay code ``C`` and its cocode ``C*``
.......................................... 

The Mathieu group ``Mat24`` operates as a permutation group on a set 
of 24 elements which we label with numbers 0,...,23 for use in Python 
and C. So it also operates on a vector space ``V = GF(2)**24``, with 
``GF(2) = {0,1}``. Here ``**`` means exponentiation. 

A vector ``v`` in a vector space over ``GF(2)`` is called a bit 
vector. We represent a bit vector as an integer, so that the 
``i``-th  bit of ``v`` (with valence ``2**i``) is the ``i``-th 
component of ``v``. 

The Golay code ``C`` ia a 12-dimensional subspace of ``V`` fixed by
``Mat24``. There are functions for checking and completing codewords 
and for getting the syndrome of a 24-bit vector.

We internally use a basis of ``V`` such that the first 12 basis
vectors are a transversal of the Golay cocode and the last 12 basis 
vectors span the Golay code. These basis vectors are listed in the
API reference.
 
We represent vectors in ``V``, ``C`` and ``C*`` and in the subset of 
octads of ``C`` as follows:

The 759 octads are numbered from 0 to 758. They do not form a vector 
space.
The 2**12 Golay code words are represented as binary numbers 0 to 4095.
The 2**12 cocode words are represented as binary numbers 0 to 4095.
A more detailed description is given in the API reference.

As usual, binary numbers representing bit vectors are added with the 
XOR operation ``^``. Unused high bits in input bit vectors are ignored. 

Functions changing an object from one representation ``xxx`` to another
representation ``yyy`` are named ``xxx_to_yyy``, where ``xxx``, 
``yyy`` is as follows::

  vect:      standard representation of a bit vector in V = GF(2)**24
             coded as a 24-bit integer. 
  vintern:   internal representation a bit vector in V as a vector in
             the basis given above, coded as 24-bit integer. 
  gcode:     representation of a Golay code word in the basis given
             given above, coded as 12-bit integer. 
  octad:     representation as an octad numbered from 0 to 758
             (in lexical order given by representation  'gcode')
  cocode:    representation as a cocode word in the basis given above, 
             coded as 12-bit integer. 


All these representations are given as unsigned integers.
 
We implement the following conversion functions::

    vect_to_vintern, vintern_to_vect, vect_to_cocode, 
    vintern_to_vect, gcode_to_vect, cocode_to_vect.

Here irrelevant bits of the input are ignored. Function
``cocode_to_vect`` returns one of many possible solutions.

In the following functions the input is checked and the function
fails in case of an error::

    vect_to_gcode, vect_to_octad, gcode_to_octad,
    octad_to_vect, octad_to_gcode

In case of failure, these C functions return a special value as 
indicated in the documentation of the function in the .c file.
The corresponding python functions raise ``ValueError`` in case
of failure. 
  
   
Function ``syndrome()`` takes a vector ``v`` and calculates its 
syndrome, which is a vector of minimum weight equivalent to ``v`` 
modulo the Golay code. Function ``cocode_syndrome()`` takes a 
``cocode`` representation of a cocode word instead.

Function ``scalar_prod()`` returns the scalar product of a Golay code
vector in ``gcode`` and a cocode vector in ``cocode`` representation.




The Mathieu group Mat24
....................... 

This class also contains support for the Mathieu group ``Mat24``.
An element of ``Mat24`` can be represented in one of the following ways::

  perm:    Representation as an array of length 24 encoding a 
           permutation of the integers 0,...,23 as a mapping.

  m24num:  Representation as an integer 0 <= i < 244823040. Here i
           is the number of the permutation in lexicographic order.
           So the identity permutation is coded as 0. 

  matrix:  Representation as a 12 x 12 bit matrix acting on the Golay
           code by right multiplication. This matrix acts on a Golay 
           code vectors (given in the 'gcode' representation) by
           right multiplication. 
           Such a matrix is implemented as an array of integers with
           each integer corresponding to a row vector of the matrix. 
           The purpose of this representation is to support 
           the Parker loop and its automorphism group. Therefore a
           row vector is implemented as a 32-bit integer.

We implement the following conversion functions::

    m24num_to_perm, perm_to_m24num, perm_to_matrix, matrix_to_perm.

There is a function ``perm_check()`` for checking if an array of 
length 24 really represents an element of the Mathieu group ``Mat24``.
All other function operating on ``Mat24`` in any way do not check if
their inputs are really in ``Mat24``. They will output garbage on bad
input, but they are not supposed to crash.

The easiest way to create a random element of ``Mat24`` is to create 
a random integer ``0 <= x < 244823040``, and to call function
``m24num_to_perm(x)``. 


Operation of the group ``Mat24`` on vectors
........................................... 

Elements of ``Mat24`` operate from the right on vectors 
``in V = GF(2)**24`` or on Golay code or cocode vectors. 
A function performing such an operation has the name::

    op_<vector>_<group>

where ``<vector>`` indicates the representation of the vector space 
and ``<group>`` indicates the representation of the group. We 
implement the functions::

    op_vect_perm, op_gcode_matrix, op_gcode_perm, op_cocode_perm.

E.g. function ``op_gcode_matrix`` operates on a Golay code word (in 
``gcode`` representation) by right multiplying an element ``m`` 
of ``Mat24`` with it. Here element ``m`` is a  12 times 12 matrix 
(in ``matrix`` representation).  


Group operation in the group ``Mat24``
......................................

Multiplication and inversion in the group ``Mat24`` is supported for 
the permutation representation ``perm``. Therefore we have functions::
 
   mul_perm, inv_perm



The Parker loop ``Pl``
......................

We support the Parker loop ``Pl`` and also its automorphism group
``AutPl``.

An element of ``Pl`` is a pair ``(v, s)``, with ``v`` a Golay code 
word and ``s`` a sign bit, as described in the API reference. We 
represent the element ``(v, s)`` as a  13-bit integer, with ``v``
given by bits 0,...,11 (in ``gcode`` representation) and the sign
``s`` given by bit 12. We call this representation of the Parker 
loop the ``ploop`` representation. So we can convert and element 
of ``C`` in 'gcode' representation to an element of ``Pl`` in 
``ploop`` representation by adjusting the sign in bit 12.

Function ``mul_ploop()`` returns the product of two elements of 
the Parker Loop. Function ``inv_ploop()`` returns the inverse of
ab element of the Parker loop.

Let ``theta`` be the cocycle for the Parker loop defined in the 
API reference. For an element ``v1`` of of ``C`` or ``Pl`` in 
``gcode`` or ``ploop`` representation, the function 
``ploop_theta(v1)`` returns the value ``theta(v1)`` (which is 
in ``C*``) in ``cocode`` representation. Function 
``ploop_cocode(v1, v2)`` returns the value of the coycle 
``theta(v1, v2)``, which is 0 or 1.


The group ``AutPl`` of standard automorphisms of the Parker loop
................................................................


An automorphism of the Parker loop is implemented as an array ``a``
of twelve 32-bit integers. The lowest 13 bits of ``a[i]`` encode 
the image of the i-th basis vector of the Parker loop. Here the 
basis of the Parker loop corresponds to the selected basis of the 
Golay code, and each basis vector has positive sign.

The bits ``13,...,24`` of the vectors ``a[i]`` encode a quadratic 
form which facilitates computations in ``AutPl``, as described 
in  section :ref:`implement-autpl-label` in the 
*Guide for developers*.

This representation of ``AutPl`` is called the ``autpl`` 
representation. We only use the ``autpl`` representation for 
elements of ``AutPl``.

Function ``perm_to_autpl(c, p)`` computes an automorphism ``m``
of the Parker loop created from an element ``p`` of ``Mat24``
(given in ``perm`` representation) and a cocode element ``c`` 
(given in ``cocode`` representation). If ``m`` is equal to
the result of ``perm_to_autpl(c, p)``, then we can get back 
``p`` and ``c`` be computing ``p = autpl_to_perm(m)`` and 
``c = autpl_to_cocode(m)``. 

Function ``cocode_to_autpl(c)`` is equivalent to function
``perm_to_autpl(c, p0)``,  where ``p0`` is the identity 
permutation. Note that::

   perm_to_autpl(c, p) = cocode_to_autpl(c) * perm_to_autpl(0, p).   

Here ``perm_to_autpl(0, p)`` is equivalent to 
*standard representative* of ``p`` in ``AutPl``, and 
``cocode_to_autpl(c)`` is a *diagonal automorphism*, as described
in section *Automorphisms of the Parker loop* of the 
API reference.

Function ``op_ploop_autpl(v, m)`` applies Parker loop automorphism 
``m`` to  element ``v`` of ``Pl`` and returns the result. 

Function ``mul_autpl(m1, m2)`` computes the product ``m1 * m2`` of 
the Parker loop automorphisms ``m1`` and ``m2``. Function 
``inv_autpl(m1)`` computes the inverse of the Parker loop 
automorphism ``m1``. 


Auxiliary functions
...................

Here is an overview of some auxiliary functions in this class. 
They are described in the corresponding function documentation:::

  bw24              bit weight of the lowest 24 bits of an integer
  lsbit24           min(24, least significant bit pos.) for an integer 
  gcode_weight      weight of a Golay code word in 'gtype' representation
  vect_to_bit_list  given a bit vector in V, it computes the lists of
                    the positions of the 0 bits and of the 1 bits of v.
  extract_b24       extract bits from bit vector using a 24-bit mask
  spread_b24        spread bit vector according to a 24-bit mask


Internal operation
..................

For switching from the standard representation to the internal
representation we use three tables with ``2**8`` entries of ``24`` bit 
length. For switching back from internal to standard representation 
we use three other tables of the same format. There are also tables 
for computing the syndrome of a vector in ``V`` with respect to the 
Golay code. There is yet another table for the cocycle ``theta`` of
the Parker loop.


Abbreviations for functions and parameters in this class
........................................................

The following list of abbreviations used in names of functions
allows to infer the action of most functions in this module::

    Abbreviation  Meaning                                   Data type

    assoc         associator (in Golay code or Parker loop)
    autpl         automorphism of the Parker loop Pl        uint32_t[12]
    bw24          bit weight of the lowest 24 bits of an int
    cap           intersection (of Golay code elements)
    cocode        element of Golay cocode C*                uint32_t
    cocycle       cocycle:  Pl times Pl  ->  {0,1}
    comm          commutator (in Golay code or Pl)
    gcode         element of Golay code C                   uint32_t
    inv           inversion (in Mat24, Pl, or AutPl)
    lsbit24       least significant bit of an integer, 
                  counting bits 0,...,23 only 
    m24num        number of an element of Mat24             uint32_t         
    matrix        element of Mat24 as binary matrix 
                  acting on the Golay code C                uint32_t[12] 
    mul           multiplication (in Mat24, Pl, or AutPl)
    net           Benes network for an element of Mat24     uint32_t[9]
    octad         number of an octad, i.e. a Golay code
                  element of weight 8                       uint32_t
    op            op_<vector>_<operation> means: 
                  apply <operation> to <vector>  
    op_all        apply operation to all vectors
    perm          element of Mat24 as a permutation         uint8_t[24]
    ploop         element of the Parker loop Pl             uint32_t    
    pow           power operator (in Pl)
    scalar_prod   scalar product (of Golay code and cocode)
    suboctad      suboctad, see function suboctad_to_cocode
    syndrome      syndrome (after decoding Golay code)      uint32_t 
    theta         cocycle theta: Pl -> C^* in Parker loop
    to            <x>_to_<y> means: return representation <y>
                  of an object given in representation <x>
    vect          vector in V = GF(2)**24                   uint32_t
    vintern       vector in V, in internal representation   uint32_t


Conventions for parameters in C functions
.........................................

Parameters of functions are either integers or arrays of integers.
Here all integer types are unsigned and of fixed length, such as
``uint8_t``, ``uint16_t`` or ``uint32_t``. 

The type of a parameter is given by a single letter in the name
of the parameter::

    Name  Meaning                                           Type
    a     array specified in documentation of function      unspecified 
    c     Golay cocode element, represented as 'cocode'     uint32_t
    m     permutation in Mat24 or automorphism of Pl
          represented as a bit matrix                       uint32_t[12]
    p     permutation in Mat24 represented as 'perm'        uint8_t[24] 
    u_<x> unsigned integer, e.g.                            unspecified
           u_exp:   integer denoting an exponent
           u_m24:   number of a permutation in Mat24
           u_octad: number of octad, 0 < u_octad < 259
           u_width: integer denoting a bit width
    v     vector in V, Golay code C or Parker loop Pl 
          represented as vect, vintern, gcode or ploop      uint32_t  

Integer input parameters have name ``u_<x>``, e.g. ``u_m24, u_exp``.
An integer computed by a function is returned as return value.
Input array parameters have a digit as a suffix, e.g.: ``v1, v2, m1``.
Output array parameters have the suffix ``_out``, e.g.: ``p_out``.
Input/output array parameters  have the suffix ``_io``, e.g.: ``m_io``.
         

Mapping C functions to python functions
.......................................


All C functions in module ``mat24_functions`` are documented.
This documentation is not repeated for the corresponding python
functions in module ``mmgroup.mat24``.

The rules for converting a C function to a python function are
as follows:

 * To obtain the name of the python function, strip off the 
   prefix ``mat24_`` from the name of a C functions.

 * To obtain the tuple of input parameters for the python function,
   take the tuple of parameters of the C function and drop are 
   parameters with suffix ``_out``. 

   In the corresponding python function, an iterable object may 
   be passed where the C function expects a pointer to an integer. 
   The minimum length of that iterable object is either clear from 
   the context or documented in the C function.

 * To obtain the return value of the python function, check if 
   the C function has any parameters with suffix ``_out``. Then 
   the sequence of returned objects is the sequence of these 
   parameters, possibly preceded by the return value if the 
   C function returns an integer value. 

   As usual in python, a sequence of length 1 is returned as a
   single object, and a sequence of length >= 1 is returned as
   a tuple. 

   A parameter with suffix ``_out`` is returned as a list of
   integers.

 * A C function may fail under certain circumstances.

   A failure of a function is indicated in the return value
   of the function. Details are given in the documentation of
   the function. The corresponding python function raises 
   ValueError if the C function fails.

 * The python function drops the return value of the C function
   from the sequence of returned python objects in certain cases.

   If the documentation of the C function contains the phrase
   'Returns 0 in case of success and (anything else) in case
   of failure' then the return value is just a status 
   indicator and hence dropped by the python function.

   If the documentation of the C function contains a phrase
   like 'Returns the length of the list ``xxx_out``' then the 
   python function adjusts the length of that returned list 
   appropriately and drops the return value.
   
 * Parameters with suffix ``_io`` refer to pointers in the C
   function and hence to iterables in the corresponding 
   python function. Here the sufix ``_io`` means that
   the function may modify that iterable object.   


"""
