.. comment

   A brief description of the Monster group


Why the Monster group is important
-------------------------------------


A **group** is a set :math:`G` equipped with a binary operation 
:math:`\cdot : G \times G \to G` (often written as multiplication, 
:math:`a \cdot b`) that satisfies the following axioms:

1. **Closure**  
   For all :math:`a, b \in G`, the product :math:`a \cdot b` is also in :math:`G`.

2. **Associativity**  
   For all :math:`a, b, c \in G` we have 

   .. math::

      (a \cdot b) \cdot c = a \cdot (b \cdot c).

3. **Identity element**  
   There exists an element :math:`e \in G` (called the *identity*) such that 
   for all :math:`a \in G` we have

   .. math::

      e \cdot a = a \cdot e = a.

4. **Inverse element**  
   For each :math:`a \in G`, there exists an element :math:`a^{-1} \in G` such that

   .. math::

      a \cdot a^{-1} = a^{-1} \cdot a = e.

If, in addition, the operation is **commutative** (:math:`a \cdot b = b \cdot a`
for all  :math:`a, b \in G`), then the group is called an **abelian group**.
In this document all groups are finite. We write :math:`|G|` for the **order**
of the group  :math:`G`, i.e. the number of elements of  :math:`G`. 

Let :math:`G` and :math:`H` be groups.  
A **group homomorphism** is a function

.. math::

   \varphi : G \to H

such that for all :math:`a, b \in G` we have

.. math::

   \varphi(a \cdot b) = \varphi(a)  \cdot  \varphi(b).
  

Here we call :math:`H` a **homomorphic image** of the group :math:`G`.
A homomorphism :math:`\varphi : G \to H` is an **isomorphism**
if it is invertible; i.e. there exist a homomorphism
:math:`\varphi^{-1} : H \to G` such that the composition of
:math:`\varphi` and :math:`\varphi^{-1}` is the identity map on
the group :math:`G`. In group theory, is usually suffices to
determine the structure of a group up to isomorphism.



A **simple group** is a group :math:`G` with :math:`|G| > 1` such that
every homomorphic image :math:`H` of :math:`G` satisfies either
:math:`|H| = |G|` or :math:`|H| = 1`. In other words, a simple group
has no nontrivial normal subgroups.

Every group can be built as a composition of simple groups. Omitting
the technical details, this construction is analogous to building a
molecule from atoms in chemistry. For a given molecule, we can specify
how many atoms of each type it contains, but such a decomposition does
not fully determine its structure. Similarly, the decomposition of a
group into simple groups provides fundamental information, though not
the complete picture.

Knowing the decomposition of a group :math:`G` into simple groups
is extremely valuable, both for analyzing the structure of :math:`G`
and for performing practical computations within it.

The classification of all finite groups (up to isomorphism) is as
infeasible as the classification of all molecules in chemistry.
However, the finite simple groups *have* been classified: by the end
of the :math:`20^{\text{th}}` century this monumental project was
completed in a series of papers spanning more than 10,000 pages.
They fall into two families:

1. **Infinite families**:

   * Cyclic groups of prime order :math:`\mathbb{Z}_p`,  
   * Alternating groups :math:`A_n` for :math:`n \geq 5`,  
   * Groups of Lie type (including classical and exceptional types over finite fields).

2. **Sporadic groups**:

   A collection of 26 exceptional groups that do not belong to the infinite families,
   the largest being the **Monster group** with about :math:`8 \times 10^{53}` elements.


While members of an infinite family of groups can often be treated uniformly,
each sporadic group requires a specialized implementation. A common approach
for computations in a sporadic group is to realize the group as a matrix
representation over a small finite field; see :cite:`Wil05` for details.

The Monster group, however, requires a representation of dimension at least
196,882, which is too large for practical computations. The second-largest
sporadic group, the Baby Monster, requires a representation of dimension 4,370,
which is feasible in practice.

The mmgroup package offers the first publicly available and computationally
efficient implementation of the Monster group, based on Conway's construction
:cite:`Con85`. While a similar approach is described in :cite:`HW03`, that
implementation is slower and not publicly accessible.


Understanding the construction of the Monster 
----------------------------------------------

The Python class |MM| supports the construction of random elements of the
Monster and computations with them. This can be done without a deep
understanding of the Monster's definition. However, our construction of the
Monster is similar to that in :cite:`Con85`. Understanding this construction
enables the user to perform considerably more sophisticated computations
within the Monster. We do not provide the full construction here, but we give
references for studying the basic structures it relies on.


1. The Mathieu group :math:`M_{24}`, the binary Golay code, and its cocode

   The Mathieu group :math:`M_{24}` is a simple group acting on
   a set :math:`\tilde{\Omega}` of 24 elements. The power set of
   :math:`\tilde{\Omega}` can be naturally interpreted as a vector space
   :math:`\mathbb{F}_2^{24}`. The group :math:`M_{24}` fixes
   a subspace of :math:`\mathbb{F}_2^{24}` of dimension 12, known
   as the *Golay code*. The Hamming distance of two different
   Golay code words is at least 8. We also consider the cocode of
   the Golay code. We adopt the notation from :cite:`CS99`, where
   all these structures are described in detail.
   

2. The Leech lattice :math:`\Lambda` and its automorphism group :math:`\mbox{Co}_0`

   The 24-dimensional Leech lattice :math:`\Lambda` is now known to be
   the lattice corresponding to the densest packing of spheres in
   dimension 24. It is described in detail in :cite:`CS99`. The
   automorphism group of :math:`\Lambda` is called :math:`\mbox{Co}_0`.
   For the  construction of the Monster, the Leech lattice modulo 2 
   (denoted by :math:`\Lambda / 2 \Lambda`) is even more important.
   The automoprhism group of :math:`\Lambda / 2 \Lambda` is a simple
   group :math:`\mbox{Co}_1`. The group  :math:`\mbox{Co}_1` is 
   obtained from :math:`\mbox{Co}_0` by identifying each matrix
   in :math:`\mbox{Co}_0` with its negative.  


3. The Parker loop  :math:`\mathcal{P}`

   After studying the Golay code, the Leech lattice, and the
   automorphism groups of these objects, the reader may consult
   Chapter 29 in :cite:`CS99`, which covers the Monster group. 
   There one will find the Parker loop as a tool for constructing the
   Monster. The Parker loop is a double cover (written multiplicatively)
   of the binary Golay code (written additively). For a deeper
   understanding of the Parker loop, see Chapter 4
   (on Symplectic 2-loops) in :cite:`Asc86`. 

   The first construction of the Monster was given by Griess :cite:`Gri82`.
   He observed that managing the signs of matrices representing generators
   of the Monster is awkward. The Parker loop provides a convenient way to
   encode these sign-related difficulties.



After studying these structures, the reader can understand Conway's
construction of the Monster :cite:`Con85` and use the Python class |MM|
to construct specific elements. Since  :cite:`Con85` does not provide
a complete set of generators of the Monster, and there are subtle
differences in sign conventions, the definitive reference for
constructing elements of the Monster is :cite:`Seysen20`.

Internally, elements of the Monster are represented as words in the
generators. Before being output, such words are reduced to a more
compact form using the highly efficient word-shortening algorithm
described in :cite:`Seysen22`.
