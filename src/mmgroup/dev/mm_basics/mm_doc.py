r"""We describe the representations of a vector in :math:`\rho_p`

The most important representation of a vector in :math:`\rho_p` is 
the *internal* representation. All operations of the monster group 
are performed on vectors in in :math:`\rho_p` in internal 
representation.

In the *internal* representation a vector is stored as an array of 
``247488`` entries, where each entry is a bit field representing 
an integer modulo ``p``. Some components of a vector are stored
twice in that array and some entries of the array are unused. This
special structure facilitates the implementation of the operations
of the monster group.

The entries of a vector are organized as tuples ``(tag, i0, i1)``. 
Here indices ``i0, i1`` refer  to a two-dimensional array as 
indicated in column ``Size`` of the  following table. As usual in 
C, entries with adjacent last index  ``i1`` are stored in adjacent 
locations. For a mathematical description of an entry 
``(tag, i0, i1)``, see section 
*The Representation of the Monster Group* in the API reference.


  .. table:: Array sizes for tags in the internal representation
    :widths: 10 24 24 24 18

    ===== ============= ============== ============= ===========
    Tag   Size          Space used     Remarks       Offset
    ===== ============= ============== ============= ===========
    ``A``   ``24 x 24``   ``24 x 32``  (1), (2)           ``0``
    ``B``   ``24 x 24``   ``24 x 32``  (1), (2), (3)    ``768``
    ``C``   ``24 x 24``   ``24 x 32``  (1), (2), (3)   ``1536``
    ``T``  ``759 x 64``  ``759 x 64``  (4),            ``2304``
    ``X`` ``2048 x 24`` ``2048 x 32``  (1),           ``50880``
    ``Z`` ``2048 x 24`` ``2048 x 32``  (1),          ``116416``
    ``Y`` ``2048 x 24`` ``2048 x 32``  (1),          ``181952``
    ===== ============= ============== ============= ===========

Remarks
  
  1. As indicated in column ``Space used``, an array of size 
     ``n`` times ``24`` is stored in a space reserved for
     an array of size ``n`` times ``32``. So the entry with
     index ``(tag, i0, i1)`` is stored at location
     ``Offset(tag) + 32 * i0 + i1``. Unused entries must
     be equal to zero.
    
  2. The entry given by ``(tag, i0, i1)`` must be equal to the 
     entry given by ``(tag, i1, i0)``. 
    
  3. The diagonal entry given by ``(tag, i0, i0)`` must be equal to 
     zero. 
   
  4. The entry with  index ``(T, i0, i1)`` is stored at location
     ``Offset(T) + 64 * i0 + i1``.

The entries of a vector in internal representation are stored
a one-dimensional array of integers of type ``uint_mmv_t``. Here
type ``uint_mmv_t`` may be one of the C integer types
``uint64_t`` or ``uint32_t``, depending on the value ``INT_BITS``.
At present ``INT_BITS`` is set to the value 64. Several adjacent
entries of a vector are stored as bit fields in a single integer
of type ``uint_mmv_t``. Entries with a lower index are stored at
bits with lower valence.


The number of bits in a bit field is always a power of two. So e.g. 
for ``p = 3`` we use ``2`` bits; for ``p = 7`` we use ``4`` bits 
with the highest bit unused. In case ``p = 2**k - 1``, legal values 
for an entry are  ``0,...,2**k - 1``, with ``2**k - 1`` equal to 
``0``. Thus negation of a value can be done by complementing all
``k`` bits of that value. Apart from negation, the matrices
corresponding to operations of the monster may add, subtract and 
half the entries of a vector (modulo ``p``). These operations can 
easily be done on several entries simultaneously by manipulating a 
just single integer of type ``uint_mmv_t``.


As indicated above, we reserve ``32`` entries for arrays of
integers modulo ``p``  with ``24`` entries. So we require about 
``25.7%`` more memory than necessary. In some cases we need this 
extra memory anyway. E.g. for ``p = 3`` a ``64``-bit integer may 
store ``32`` entries, so that there will always be a slack of 
``8`` entries when storing ``24`` entries. 

Function ``mm_aux_mmv_size(p)`` returns the number of integers 
of type  ``uint_mmv_t`` required for storing a vector in external
representation. 

When writing or calculating entries with tags ``A, B, C`` then
the high-level function for manipulating vectors in internal 
representation make sure that e.g. entries with indices
``(A, i0, i1)`` and ``(A, i1, i0)`` are always set to the same 
value. These functions also make sure that unused bits in a bit 
field and unused bit fields are set to zero. Note that a bit field 
with value zero may contain the value ``p`` instead.



The *external* representation of a vector in :math:`\rho_p`

There is also a so-called *external representation* of a vector in R_p. 
This is used to facilitate the access to vectors by external modules. 
Here the vector is represented as an array of 196884 integers of type
uint8_t. Basis vectors are ordered similar to the ordering for the
internal representation, but here the entries are in one-to-one
correspondence with the basis vectors. In the external representation 
there are no unused or duplicated entries.

More precisely, the order of the entries is:

  .. table:: Order of entries in external representation
    :widths: 25 25 25 25

    ===============  ===========  ==============  ===========
    Entries          Condition    No of entries   Offset
    ===============  ===========  ==============  ===========
    ``(A, i0, i1)``  ``i0 = i1``        ``24``         ``0``
    ``(A, i0, i1)``  ``i0 > i1``        ``276``       ``24``
    ``(B, i0, i1)``  ``i0 > i1``        ``276``      ``300``
    ``(C, i0, i1)``  ``i0 > i1``        ``276``      ``576``
    ``(T, i0, i1)``                  ``759*64``      ``852``
    ``(X, i0, i1)``                 ``2048*24``    ``49428``
    ``(Z, i0, i1)``                 ``2048*24``    ``98580``
    ``(Y, i0, i1)``                 ``2048*24``   ``147732``
    ===============  ===========  ==============  ===========

Indices (``tag, i,j``) for ``tag = A, B, C``, ``i > j`` are ordered 
as follows:

   *  ``(1,0),``
   *  ``(2,0), (2,1),``
   *  ``(3,0), (3,1), (3,2),``
   *  ``...``
   *  ``(i,0), (i,1), ..., (i,i-1),``
   *  ``...``
   *  ``(24,0), (24,1), ..., (24,23).``

Function ``mm_aux_bytes_to_mmv()`` converts a vector from external  
to internal representation, Function ``mm_aux_mmv_to_bytes()`` does 
the inverse conversion.



The *sparse* representation of a vector in :math:`\rho_p`

The Python interface to vectors in :math:`\rho_p` is optimized
for readability and not for speed. Here a typical task is to read
and modify single entries of a vector. In the internal 
representation the coordinates of a vector are indexed by tuples
containing a string. Transferring tuples or strings from Python
to C is awful. Here we need a representation of a vector in
:math:`\rho_p` where a single entry of a vector can be stored
in an integer variable.

A vector in :math:`\rho_p` can be stored in the *sparse* 
representation. 
Here a vector is stored as an array of 32-bit integers, where each 
entry stands for a multiple of a basis vector. A component of
a vector is stored in the bit fields of an integer as a tuple
``(tag, i0, i1, value)``. Here the tuple  ``(tag, i0, i1)`` is as
in the external representation, and ``value`` is the value of the
coordinate of the vector corresponding to ``(tag, i0, i1)``.
Entries with coordinate zero may be dropped. A 32-bit integer
encodes a tuple ``(tag, i0, i1, value)`` in bit fields as shown 
in the following table.


   .. table:: Bit fields in the sparse representation
     :widths: 20 80

     ======   ======================================================
     Bits     Meaning
     ======   ======================================================
     27..25   ``Tag``: 
              ``A = 1``, ``B = 2``, ``C = 3``, ``T = 4``,
              ``X = 5``, ``Z = 6``, ``Y = 7``
     24..15   Index ``i0``
     13.. 8   Index ``i1``
      7.. 0   The ``value`` of the coordinate of the basis vector;  
              if the modulus ``p`` is ``2**k - 1``  then only 
              the lowest ``k`` bits are evaluated.
     ======   ======================================================
            
In a C function the length of a sparse representation of a vector 
must be given as a parameter to the function. 
The order of the entries is irrelevant in the sparse 
representation. A sparse representation generated by a C function
contains at most one entry for each tuple ``(tag, i0, i1)``.
On input, entries with several equal tuples ``(tag, i0, i1)`` are 
accepted. Unless stated otherwise, the corresponding values of
such equal tuples are added.

In an entry with tag ``A``, ``B``, or ``C`` generated by this module
we always have ``i0 >= i1``. The ``value`` of an entry generated by 
this module is always less than the modulus ``p``.

When reading an entry, a ``value`` with ``0 <= value <= p`` is
accepted. Entries with tag ``A``, ``B``, or ``C`` and ``i < j`` are 
also accepted. Illegal tags or indices are usually ignored on input.

"""