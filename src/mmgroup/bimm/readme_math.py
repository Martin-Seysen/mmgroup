r"""The mapping from the Coxeter group into the Bimonster. 
-----------------------------------------------------------

In this section we briefly describe the mathematical tools used
for the mapping from the Coxeter group :math:`\mbox{IncP3}` to the 
Bimonster.

Therefore we use the Norton's presentation :cite:`Nor02` of the 
Monster. Then we follow the ideas in Farooq's thesis :cite:`Far12` 
for extending that presentation of the Monster to a homomorphism from 
the Coxeter group :math:`\mbox{IncP3}` to the Bimonster.
We assume that the reader is familiar with :cite:`Nor02` and 
we will adopt notation from ibid.


Norton's presentation of the Monster and the Bimonster
.......................................................

Let :math:`\mbox{IncP3}` be the Coxeter group generated by the 
nodes of the projective plane :math:`\mbox{P3}` as above, and let
:math:`\mbox{AutP3}` be the automorphism group of :math:`\mbox{P3}`.
Let :math:`\,\mbox{IncP3}:\mbox{AutP3}\,` be the extension group 
where :math:`\mbox{AutP3}`  acts naturally on the generating
reflections of the Coxeter group :math:`\mbox{IncP3}`.

Norton :cite:`Nor02` defines a presentation :math:`(s,t,u,v,x,\alpha)`
of the group :math:`\mbox{IncP3}:\mbox{AutP3}`, where 
:math:`s,t,u \in \mbox{AutP3}`, :math:`v = \prod_{i=1}^{12} P_i`,  
:math:`x = P_0 L_0`, :math:`\alpha = P_0 v`. 
Then he adds a relation that maps :math:`\mbox{IncP3}` to the 
Bimonster :math:`\mathbb{M} \wr 2`, thus obtaining a presentation of
:math:`(\mathbb{M} \wr 2):\mbox{AutP3}`. Next he shows that
this is actually a presentation of the direct product 
:math:`(\mathbb{M} \wr 2) \times \mbox{AutP3}`. Then he adds a
set of relations that cancel the factor :math:`\mbox{AutP3}`, 
leading to a presentation of the Bimonster with the generators
:math:`(s,t,u,v,x,\alpha)`. It turns out that all these generators,
except for :math:`\alpha`, are in the subgroup
:math:`(\mathbb{M} \times \mathbb{M})` of index 2 of the Bimonster
:math:`\mathbb{M} \wr 2`. Finally, he gives a further set of 
relations in the generators :math:`(s,t,u,v,x)` cancelling the second 
factor of the direct product :math:`(\mathbb{M} \times \mathbb{M})`,
thus obtaining a presentation of the Monster :math:`\mathbb{M}`.



TODO: Add more documentation!!!
"""

