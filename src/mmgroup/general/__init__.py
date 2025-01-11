r"""We support abtract groups acting on certain sets.

.. warning::

  The functions in this module are yet experimental and subject to
  change in the future. We recommend to contact the author  before
  using them in research projects.

Let :math:`G` be a finite group and :math:`V` be a set. An action
of :math:`G` on :math:`V` is a mapping :math:`\rho` from :math:`G`
into a permutation group of the set :math:`V`. The current version
supports actions on sets :math:`V` with a certain structure only.
Here :math:`V` may be

- A vector space :math:`V = \mbox{GF}_2^n, n \leq 24`, with  :math:`G`
  acting as a group of linear transformations on :math:`V`.
  Such an action is modelled by class ``Orbit_Lin2``.

- An elementary Abelian 2 group :math:`V` of structure
  :math:`2^n, n \leq 32`, with :math:`G` acting on :math:`V` via a
  homomorphism from :math:`G` to :math:`V`.
  Such an action is modelled by class ``Orbit_Elem2``.

In both cases a group is given by a set of generators
:math:`(g_{0} \ldots g_{k-1})` and a mapping :math:`\rho` that
maps each of these generators to an element of :math:`V`. The only
requirements for the elements :math:`g_{i}` are that they support
multiplication (for the group operation) and exponentiation with
an integer exponent (e.g. for computing inverses an the neutral
element). The mapping :math:`\rho` must be implemented as a
function taking a single group element :math:`g_{i}` and returning
:math:`\rho(g_{i})` as an element of :math:`V`. The encoding of
elements of :math:`V` is is described in the documentation of the
corresponding class.

The methods of these classes allow the computation of orbits of
:math:`V` under the action of :math:`G`. We also may compute
Schreier vectors on an orbits of :math:`V` (under the action of
:math:`G`), as described in :cite:`HE05`, Section 4.1. If :math:`G`
is acting on an elementary Abelian 2 group :math:`V`, we may
compute the kernel of that action instead.

The two types of actions mentioned above are desigend to
compute in the managable 2-local subgroups of the Monster, especially
in the groups :math:`G_{x0}` of structure :math:`2^{1+24}.\mbox{Co}_1`,
and in the group :math:`N_{0}` of structure
:math:`2^{2+11+22}.(\mbox{M}_{24} \times \mbox{S}_{3})`.

In principle the functionality in this module can be used for
computing the order of a subgroup of such a managable subgroup,
or to perform a constructive membership test for such a subgroup.
The current implementation should be considered as
*under construction* and may be extended in future versions.

Class ``Random_Subgroup`` in this module implements a simple
random generator for generating pseudorandom elements of a
group given by generators.
"""

from mmgroup.general.orbit_lin2 import Orbit_Lin2
from mmgroup.general.orbit_lin2 import Orbit_Elem2
from mmgroup.general.orbit_lin2 import Random_Subgroup
