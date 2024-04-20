

.. _subgroup-mat24-label:

Subgroups of the Mathieu group :math:`M_{24}` 
---------------------------------------------


An amalgam of many (mostly 2-local) large subgroups of the  Monster
:math:`\mathbb{M}` has been constructed in :cite:`Iva09`. For each
of these subgroups (and, hopefully, also for their intersections) 
we want an algorithm that constructs an (almost) uniform random
element of that group. This is a complicated task, which might not
even be possible for all cases. As a first step towards this goal 
we will describe the intersections of some subgroups of the Mathieu 
group :math:`M_{24}`. These subgroups of :math:`M_{24}` are involved
in the subgroups of :math:`\mathbb{M}` in the amalgam.
We will also give some information how to construct ramdom elements
of these subgroups of :math:`M_{24}` and of their intersections.
 
TODO: This section is yet a fragment and will be completed in a 
future version. 


Some subgroups of the the Mathieu group :math:`M_{24}` 
......................................................

We define some subgroups of the Mathieu group :math:`M_{24}`
as the stabilizers of certain subsets (or sets of subsets) of the 
set :math:`\tilde{\Omega} = \{0,\ldots,23\}`, as shown in the 
following table:



.. math::
    \begin{array}{|c|c|c|c|}
    \hline 
    \mbox{Name} & \mbox{Mnemonic} & \mbox{Stabilizes the set} &
       \mbox{Structure}  \\
    \hline
    M_2 & \mbox{2-set} & \{2,3\} &
       M_{22}:2 \\ 
    M_o & \mbox{octad}  &  \{0,\ldots,7\} &
       2^4:A_8  \\
    M_t & \mbox{trio}   & 
       {\{ \{8i,\ldots,8i+7\} \mid i < 3 \}} &
       2^6:(S_3 \times L_3(2))   \\
    M_s & \mbox{sextet} & 
       {\{ \{4i,\ldots,4i+3\} \mid i < 6 \}} &
       2^6:3.S_6 \\
    M_l & \mbox{(line)} & \{ \{2i, 2i+1\} \mid 4 \leq i < 12 \} &
       2^{1+6}:L_3(2)  \\
    M_3 & \mbox{3-set} & \{1, 2,3\}   &
       L_3(4):S_3 \\
    \hline 
    \end{array}

All these groups, except for :math:`M_l`, are maximal subgrous of
:math:`M_{24}` and discussed in :cite:`CS99`, Ch. 11. The structure
of :math:`M_l` (as a subgoup of :math:`M_o`) is described in
:cite:`Iva09` in the text after Lemma 4.1.3. The central element of
:math:`M_l` of order 2 is obtained by exchanging all entries
:math:`2i` with :math:`2i+1` for  :math:`i \geq 4`. We have the
following inclusions:

.. math::
    \begin{array}{lll}
   M_l \subset M_o \, , \; &
   M_l \cap M_2 \subset M_t \, , \; &
   M_t \cap M_2 \subset M_o \cap M_l \, , \\
   M_3 \cap M_l \subset M_s \, , \; &  
   M_3 \cap M_t \subset M_o \cap M_s \, , \; &  
   M_3 \cap M_t \cap M_l \subset M_2 \; . \\
   \end{array}

Membership of an element of :math:`M_{24}` in any of these subgroups
can easily be checked computationally. Thus an inclusion of two 
intersections of these subgroups can be disproved computationally by 
finding an element of :math:`M_{24}` contadicting that inclusion. 
This computation is done in function ``test_mat24_rand`` in module 
``mmgroup.tests.test_mat24`` during the standard test. 


In the remainder of this subsection we will prove the inclusions
given above. The first inclusion is obvious. Next we justify the 
following two inclusions. 


By :cite:`CS99`, Ch. 11, the group :math:`M_o` acts on the
octad :math:`o = \{0,\ldots,7\}` as the alternating group
:math:`A_8`, and on the complement :math:`\bar{o}` of :math:`o`
as the affine group on :math:`\mathbb{F}_2^4`. Here the affine
structure of :math:`\bar{o}` is given by the last 4 digits of
the binary representations of its entries. Fixing e.g. the
element 8 of :math:`\bar{o}`, the set  :math:`\bar{o}` acquires 
the structure of the linear space :math:`\mathbb{F}_2^4`; and
we obtain the well-known isomorphism :math:`L_4(2) \cong A_8`.
Subspaces of dimension 1, 2, and 3 of the linear 
space :math:`\mathbb{F}_2^4` will be called lines, planes,
and hyperplanes.

We have following correspondences between  :math:`o` and its
complement :math:`\bar{o}` 


.. math::
    \begin{array}{|l|l|l|}
    \hline
    \mbox{Linear space } \bar{o} =  \mathbb{F}_2^4  
       & \mbox{Octad } o  &
          \mbox{Remarks} \\
    \hline
    \mbox{line} & \mbox{Steiner system }  S(3,4,8)  &
        \mbox{defines an affine structure on } o \\
    \hline
    \mbox{plane} & \mbox{sextet refining } o  &
        \mbox{tetrads of the sextet in } \bar{o} \mbox{ are} \\
        & & \mbox{planes parallel to that plane} \\
    \hline
    \mbox{hyperplane} & \mbox{Steiner system }  S(3,4,8)  &
        \mbox{defines an affine structure on } o 
    \mbox{} \\
    \hline
    \mbox{symplectic form} & \mbox{set of 2 elements}   &
         \\
    \hline
    \mbox{line incident} & \mbox{partition into 4} & \\
        \mbox{with hyperplane}  & \mbox{sets of 2 elements}  \\
    \hline
    \end{array}


These correspondences (except for the last) are stated in
:cite:`Atlas`. For the last correspondence we refer to 
:cite:`CS99`, Ch. 10.1.4. The tetrads of a Steiner system
corresponding to a line or a hyperplane are the tetrads in the 
sextets corresponding to the planes incident with that line
or hyperplane; here we take tetrads that are subsets of 
octad :math:`o` only. Methods for computing with these objects 
are presented in :cite:`CS99`, Ch. 11, and implemented in file 
``mat24_functions.c``.

Let :math:`f` be the symplectic form on :math:`\mathbb{F}_2^4`
corresponding to a pair :math:`y` of elements of :math:`o`.
For two different nonzero elements :math:`x_1, x_2` of 
:math:`\mathbb{F}_2^4` let :math:`s` be the sextet corresponding
to the plane generated by :math:`x_1`  and :math:`x_2`, and
let :math:`\tau` be a tetrad in  :math:`s`  which is a subset
of  :math:`o` . Then the scalar product  :math:`x_1`  and 
:math:`x_2` with respect to :math:`f` is the size of the set
:math:`y \cap \tau` (modulo 2).
 
Let :math:`l' = \{8,9\}, s'=  \{8,9,10,11\}` and
:math:`t' = \{8,\ldots,15\}`. Then the line  :math:`l'`, the
plane   :math:`s'`, and the hyperplane  :math:`t'` are
mutually incident. A staightforward calculation using the
facts stated above shows
that  :math:`t'` is the orthogonal complement of  :math:`l'`
with respect to the symplectic form  :math:`f'` corresponding
to the subset :math:`\{2,3\}` of :math:`o`. 
Obviously, :math:`M_t \cap M_2 subset M_o`. So this proves
:math:`M_l \cap M_2 \subset M_t` and
:math:`M_t \cap M_2 \subset M_o \cap M_l`.


Next we will show the last three inclusions listed above.
The group :math:`M_3` fixes the set :math:`y= \{1,2,3\}`.
We have  :math:`M_3 \cap M_l \subset M_o`, so 
:math:`M_3 \cap M_l` fixes a Steiner system in :math:`o`,
an it turns out the the set :math:`\{0,1,2,3\} \supset y` 
is in that Steiner system. This proves 
:math:`M_3 \cap M_l \subset M_s`. We obviously have
:math:`M_3 \cap M_t \subset M_o`, and hence
:math:`M_3 \cap M_t` fixes another Steiner system in :math:`o`.
Here the set :math:`\{0,1,2,3\}` is also in this Steiner system,
implying :math:`M_3 \cap M_t \subset M_o \cap M_s`. 

The group  :math:`M_t \cap M_l` fixes a line incident with a 
hyperplane in :math:`\mathbb{F}_2^4`. Hence it fixes a partition
of :math:`o` into four unordered pairs, and it turns out that
the set :math:`u` of these pairs is  
:math:`\{ \{2i, 2i+1\} \mid 0 \leq i < 4 \}`.
Thus  :math:`M_3 \cap M_t \cap M_l` fixes :math:`u` and
:math:`y`, and hence also the set :math:`\{2,3\}`. This proves
:math:`M_3 \cap M_t \cap M_l \subset M_2`.


Some 2-local subgroups of the Monster
.....................................

The subgroups of  :math:`M_{24}` constructed above are involved
in some large subgrups of the Monster as described in the
following table. 


.. math::
    \begin{array}{|c|c|c|c|}
    \hline 
    \mbox{Name} & \mbox{Involved in subgroup of} \,  \mathbb{M} \\
    \hline
    M_2 &    H^+ = 2.\mbox{B} \\ 
    M_o &   G_{10}  = 2^{10+16}.O_{10}^+(2)  \\
    M_t &   G_5^{(t)}  = 2^{5+10+20}.(S_3 \times L_5(2))  \\
    M_s &  G_3 = 2^{3+6+12+18}.(L_3(2) \times L_3(2))   \\
 %   M_l &  G_5^{(l)}   \\
    M_3 &  2^2.\!\hphantom{|}^2\!E_6(2):S_3 \\
    \hline 
    \end{array}




TODO: More details will be given in a future version.


