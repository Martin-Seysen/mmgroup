.. comment

   This ist the content of the following chapter in the API refernece:
   Generating random elements of the Monster and of some of its subgroups

In this section we present various subgroups of the Monster known by 
the ``mmgroup`` package. Each of these groups is labelled ba a string.
In the documentation of class |MM| we have already seen the basic 
way how to generate a random element of a subgroup of the Monster.
E.g. calling method ``MM('r', 'G_x0')`` generates an element of the
subgroup :math:`G_{x0}` of structure  :math:`2^{1+24}.\mbox{Co}_1`.
We may use some other names instead of ``'G_x0'`` for generating
random elements in other subgroups of the Monster. Details will be
given in the following subsections.

At present, these names are used to generate random elements of
the corresponding subgroup. In future versions of ``mmgroup``
we may also deal with constructive membership testing in these 
groups, and with intersections of these groups. So we will only
label groups, where the embedding of this group in the Monster is
sufficiently well understood.  



Describing a random subgroup of the Monster
.............................................


The subgroups of the Monster recognized by the ``mmgroup`` package
fall in two (not necessarily disjoint) categories. There are small 
subgroups that are given by the sets of their generators; and there 
are large 2-local subgroups which are given by their centralizers. 
Each subgroup has one or more names; and it is described in the 
following two sections.  

The intersection of two groups is obtained by joining their names
with an ``'&'`` character. E.g, ``'G_x0 & N_0'`` means the
intersection of the two groups with names  ``'G_x0'`` and ``'N_0'``.
We have chosen the definition of the subgroups carefully, so that 
at least in principle their intersections can be computed. 

But not all possible computations have been implemented. For each
group there is a column 'Tested' in the following tables, with
an entry 'full', 'yes', or 'no'. If all groups in a tuple are
supported 'full' then we may compute intersections of these groups.
If a group has an entry 'yes' then we may use that group. Otherwise
using that group might lead to an error. 

A function dealing with the name of a subgroup  may raise 
``NotImplementedError`` if that subgroup is not supported. 


Caution:

The user should check the corresponding entry 'Tested' in the tables
in the follwing two subsections before using the name of a group!
Using a name of a group that is not supported might fail without notice!
  



Some small subgroups of the Monster
.............................................


The following small subgroups of the Monster are recognized:

.. table:: Some small subgroups of the Monster
    :widths: 15 30 40 15 

    ============== ============================================ ========================================================================== ======
    Name           Structure                                    Generators                                                                 Tested
    ============== ============================================ ========================================================================== ======
    N_0            :math:`2^{2+11+22}.(S_3 \times M_{24})`      :math:`x_d, y_d, x_\delta, x_\pi, \tau`                                      full
    N_0_e          :math:`2^{2+11+22}.(3 \times M_{24})`        :math:`x_d, y_d, x_\delta, x_\pi, \tau, \delta \;  \mbox{even}`              full 
    N_x0           :math:`2^{1+24+11}.M_{24}`                   :math:`x_d, y_d, x_\delta, x_\pi`;  equal to :math:`G_{x0} \cap N_0`         full
    N_x0_e         :math:`2^{2+11+22}.M_{24}`                   :math:`x_d, y_d, x_\delta, x_\pi, \delta \; \mbox{even}`                     full
    Q_x0           :math:`2^{1+24}`                             :math:`x_d, x_\delta`                                                        no
    AutPL          :math:`2^{11+1}.M_{24}`                      :math:`x_\delta, x_\pi`                                                      no
    AutPL_e        :math:`2^{11}.M_{24}`                        :math:`x_\delta, x_\pi, \delta \; \mbox{even}`                               no
    ============== ============================================ ========================================================================== ======


Notation for the generators is as in  :py:mod:`~mmgroup.mm_group`.
 


Large 2-local subgroups of the Monster
.............................................

We have implemented (or will implement) a number of 2-local subgroups
of the Monster. Each of these has a name *G* as shown in the following
table. The group is defined as the normalizer of an elementary Abelian
2 group  *Z* = *Z(G)*. In the following table we will list one or
more names for such a group, and the structure of that group, if known.
We also display a generating system for the group *Z(G)* and the
exponent ``e`` of the order :math:`2^e`of  *Z(G)* . The basic reference
for the constrution of these 2-local subgroups is :cite:`Iva09`; but we
also support the names of thise groups geven in  :cite:`Con85`
 
The implementation of these large 2-local subgroup is yet experimental!  



.. table:: Large 2-local subgroups of the Monster
    :widths: 18 29 33 4 11 

    ============== ============================================ ========================================================================== === ======
    Name           Structure                                    Generators of the 2 group Z                                                  e Tested
    ============== ============================================ ========================================================================== === ======
    G_1, G_x0      :math:`2^{1+24}.\mbox{Co}_1`                 :math:`x_{-1}`                                                               1  full
    G_2, N_0       :math:`2^{2+11+22}.(S_3 \times M_{24})`      :math:`x_{-1}, x_{\pm \Omega}`                                               2  full
    G_3            :math:`2^{3+6+12+18}\!.(L_3(2) \times 3S_6)` :math:`x_{-1}, x_{\pm \Omega}, x_\omega`                                     3   no
    G_4            ?                                            :math:`x_{-1}, x_{\pm \Omega}, x_\omega, x_{\{0,1,4,5\}}`                    4   no
    G_5t           :math:`2^{5+10+20}.(S_3 \times L_5(2))`      :math:`x_{-1}, x_{\pm \Omega}, x_\omega, x_{\{0,1,4,5\}}, x_{\{0,2,4,6\}}`   5   no
    G_5l           ?                                            :math:`x_{-1}, x_{\pm \Omega}, x_\omega, x_{\{0,1,4,5\}}, x_{\{0,2,4,7\}}`   5   no
    G_10           :math:`2^{10+16}.\Omega_{10}^+(2)`           :math:`x_d, y_d, x_{ij}, \quad 0 \leq i < j < 8,`  
                                                                :math:`d \in \{-1, \pm \Omega, d_8\}`                                       10   no
    B              :math:`2.B`                                  :math:`x_{\{2,3\}}`                                                          1  yes
    2E_6           :math:`2^2.{}^2E_6(2):S_3`                   :math:`x_{\{2,3\}}, x_{\{1,2\}}`                                             2   no
    ============== ============================================ ========================================================================== === ======


Notation is as in the previous section. Furthermore,
:math:`d_8` is (a preimage in the Parker loop of)
the standard octad :math:`\{0,1,2,3,4,5,6,7\}`.
:math:`\omega` is the element of the Golay cocode given by 
the tetrad :math:`\{0,1,2,3\}` as in :cite:`Seysen20`. 


The names *G_x0* and *N_0* are adopted from :cite:`Con85`. The names and
structures of *G_1, G_2, G_3, G_5t, G_5l* and *G_10* are taken from
:cite:`Iva09`, Chapter 4.1 et seq. We have *G_1 = G_x0* and
*G_2 = N_0*. The other names indicate the large simple subgroups
involved in that group. Note that all these subgroups are
maximal in the Monster, except for *G_4* and *G_5l*. 


Our chosen generating system of the Monster contains two generators
:math:`\tau` and :math:`\xi`, and a set :math:`\Gamma_1` of some
more generators. The subgroup  *G_1* is generated by :math:`\xi` and
:math:`\Gamma_1`; and the subgroup *G_2* is generated by :math:`\tau`
and :math:`\Gamma_1`. Most of the difficulty in computing in the
monster group comes from the interaction of the generators :math:`\xi`
and :math:`\tau`. So computing in *G_1*  or  *G_2* is much easier than
computing in the Monster.
 

Our large 2-local subgroups *G* have been chosen so that all of them
(except for *G_1* and *G_2*) contain both generators
:math:`\tau` and :math:`\xi`. The intersection of such a group *G*
with *G_2* (of structure :math:`2^{2+11+22}.(S_3 \times M_{24})`)
has a subgroup :math:`H` of :math:`M_{24}` as a factor group. These
factor groups have been deacribed in considerable detail in
:cite:`Iva09`,  Chapter 4.1 and 7. Such a subgroup :math:`H` of
:math:`M_{24}` can be described by its action on the set
:math:`\Omega` on which :math:`M_{24}` acts naturally. In the
following table we characterise the groups  :math:`H` by the
subsets (or partitions) of :math:`\Omega` fixed by the group. 



.. table:: Subgroups :math:`M_{24}` involved in the large 2-local subgroups
    :widths: 10 25 65

    ====== ==================================== ==========================================================================================================
    Name   Subgroup :math:`H` of :math:`M_{24}` Structure in set  :math:`\Omega` fixed by  :math:`H`  
    ====== ==================================== ==========================================================================================================
    G_3    :math:`2^{6}: 3.S_6`                 sextet  :math:`\{\{4i,\ldots,4i+3\} \mid 0 \leq i < 6 \}` 
    G_4    ?                                    partition    :math:`\{\{2i, 2i+1\} \mid  0 \leq i < 12 \}`  
    G_5t   :math:`2^{6}.(S_3 \times L_3(2))`    trio   :math:`\{\{8i,\ldots,8i+7\} \mid 0 \leq i < 3 \}`
    G_5l   ?                                    octad :math:`d_8`, partition :math:`\{\{2i, 2i+1\} \mid  4 \leq i < 12 \}` of :math:`\Omega \setminus d_8`
    G_10   :math:`2^{4}.A_8`                    octad :math:`d_8 = \{0,1,2,3,4,5,6,7\}`                                    
    B      :math:`M_{22}:2`                     duad :math:`\{2,3\}`  
    2E_6   :math:`L_3(4):S_3`                   triad :math:`\{1,2,3\}` 
    ====== ==================================== ==========================================================================================================



