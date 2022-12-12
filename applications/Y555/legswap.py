import inc_p3
from bimm import BiMM, AutP3, AutP3_BiMM, P3_node
from mmgroup import MM


legswap_autp3 = AutP3('a:a, c1:c1, c2:c3, c3:c2')
legswap = AutP3_BiMM(legswap_autp3)
m1, m2, e = legswap.decompose()
assert m1 == m2
print(m1)
I1, _ =  m1.conjugate_involution()
I2, _ =  m2.conjugate_involution()
print(I1, I2)