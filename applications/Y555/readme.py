r"""Find a homomorphism from the Coxeter group :math:`Y_{555}` into the Bimonster


The author has coded a (yet undocumented) implementation of this
homomorphism on his local computer. This will be ported to the 
github repository in a future version.

Warning: The following ducumentation is just as stub!


Description of the task to be done
..................................


The wreath product :math:`\mathbb{M} \wr 2` is of the Monster
:math:`\mathbb{M}` with the group :math:`Z_2` of order 2 is called 
the *Bimonster*. In ATLAS notation :cite:`Atlas` the Bimonster has 
structure :math:`(\mathbb{M} \times \mathbb{M}).2`. It is generated
by :math:`\mathbb{M} \times \mathbb{M}` and an involution
swapping the two copies of  :math:`\mathbb{M}`
 
 

A Coxeter group is a group generated by a finite set 
:math:`s_1,\ldots, s_k` of involutions with a set of relations
:math:`(s_i, s_j)^{\alpha_{i,j}} = 1` for each pair :math:`s_i, s_j`
of involutions. The Coxeter groups discussed in this application
are given by graphs, where the vertices of the graph correspond 
to the generators :math:`s_i`. If two vertices  :math:`s_i, s_j` 
are connected by an edge then we have a relation
:math:`(s_i, s_j)^3 = 1`; otherwise we have a relation 
:math:`(s_i, s_j)^2 = 1`. In the last case the generators
:math:`s_i` and  :math:`s_j` commute. 


Norton and Ivanov have shown that the Bimonster is isomorphic
to a certain Coxeter group, together with a single 
additional relation, see :cite:`Nor92`, :cite:`Iva92`.
That presentation of the Bimonster is called :math:`Y_{555}`.

The graph corresponding to the Coxeter relations in
:math:`Y_{555}` is given in the following Figure 1. The 
names of the generating refelctions of :math:`Y_{555}`
in that figure are as in  :cite:`Atlas`.


.. tikz:: Coxeter relations in the group Y_555

   [
   dot/.style = {circle, fill, minimum size=#1,
   inner sep=0pt, outer sep=0pt},
   dot/.default = 4pt  % size of the circle diameter 
   ]  
   \node[dot](f1) at (0,5){} ;
   \node [left=1mm of f1] {$f_1$};
   \node[dot](e1) at (0.5,4.5){};
   \node [left=1mm of e1]{$e_1$};
   \node[dot](d1) at (1,4){};
   \node [left=1mm of d1]{$d_1$};
   \node[dot](c1) at (1.5,3.5){};
   \node [left=1mm of c1]{$c_1$};
   \node[dot](b1) at (2,3){};
   \node [left=1mm of b1]{$b_1$};
   \node[dot](a) at (2.5,2.5){};
   \node [above=1mm of a]{$a$};
   \path[-]
   (f1) edge (e1) edge (d1) edge (c1) edge (b1) edge (a)	  ;
   \node[dot](b3) at (2.5,2.0){} ;
   \node [left=1mm of b3] {$b_3$};
   \node[dot](c3) at (2.5,1.5){} ;
   \node [left=1mm of c3] {$c_3$};
   \node[dot](d3) at (2.5,1.0){} ;
   \node [left=1mm of d3] {$d_3$};
   \node[dot](e3) at (2.5,0.5){} ;
   \node [left=1mm of e3] {$e_3$};
   \node[dot](f3) at (2.5,0.0){} ;
   \node [left=1mm of f3] {$f_3$};
   \path[-]
   (f3) edge (e3) edge (d3) edge (c3) edge (b3) edge (a)	  ;
   \node[dot](b2) at (3,3){} ;
   \node [right=1mm of b2] {$b_2$};
   \node[dot](c2) at (3.5,3.5){} ;
   \node [right=1mm of c2] {$c_2$};
   \node[dot](d2) at (4,4){} ;
   \node [right=1mm of d2] {$d_2$};
   \node[dot](e2) at (4.5,4.5){} ;
   \node [right=1mm of e2] {$e_2$};
   \node[dot](f2) at (5,5){} ;
   \node [right=1mm of f2] {$f_2$};
   \path[-]
   (f2) edge (e2) edge (d2) edge (c2) edge (b2) edge (a)	  ;	  


We let :math:`Y_{555}` be the Coxeter group given by the generators
and relations in Figure 1 together with the following additional 
relation:

.. math:: 

   (a b_1 c_1 a b_2 c_2 a b_3 c_3)^{10} = 1 \, . 


The purpose of this application is to implement the mapping
from the group  :math:`Y_{555}` to the Bimonster. In the sequel
we write  :math:`Y_{555}` for both, the graph in Figure 1, and
the group defined above.


Extending :math:`Y_{555}` to the incicdence graph of a projective plane
----------------------------------------------------------------------- 


The graph :math:`Y_{555}` can be extended to the incidence graph 
:math:`\mbox{IncP3}` of the projective plane :math:`\mbox{P3}` 
over the field :math:`\mathbb{F}_3`. The 26-Node Theorem states 
that there is a homomorphism from the Coxeter group given by
:math:`\mbox{IncP3}` to the Bimonster :math:`\mathbb{M} \wr 2`,
see :cite:`CNS88`. Defining relations in the generators 
of that Coxeter group defining that homomorphism are well
known see e.g. :cite:`Nor02`, :cite:`Far12` for an 
overview. We also write :math:`\mbox{IncP3}` for the Coxeter
group corresponding to  to the incidence 
graph :math:`\mbox{IncP3}`.

In practice it turns out that working with  :math:`\mbox{IncP3}` 
is more flexible than working with :math:`Y_{555}`. 

The projective plane  :math:`P3` contains 13 points :math:`P_i` and
13 lines  :math:`L_i`, :math:`0 \leq i < 13`. Point  :math:`P_i` is
incident with line :math:`P_j` if :math:`i + j` is 0, 1, 3, or 9
(modulo 13). This construction of :math:`P3` is also used in
:cite:`CNS88`,  :cite:`Nor02`, and :cite:`Far12`. 
According to :cite:`Far12`, we may map the 16 vertices
of :math:`Y_{555}` to a subset
of the set of 26 vertices of :math:`\mbox{IncP3}` as follows:


.. math::

    \begin{array}{c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c}
    a & b_1 & b_2 & b_3 & c_1 & c_2 & c_3 & d_1 & d_2 & d_3 & 
      e_1 & e_2 & e_3 & f_1 & f_2 & f_3  \\
    \hline
    P_0 & L_3 & L_0 & L_1 & P_{10} & P_1 & P_2 & 
    L_4 & L_2 & L_7 & P_5 & P_7 & P_6 &
    L_8 & L_6 & L_{10}
    \end{array}


Implementing the Bimonster and the homomorphism from :math:`IncP3` to it
........................................................................


TODO: Yet to be described!!!

"""

