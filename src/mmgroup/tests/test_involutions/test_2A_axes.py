from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint #, shuffle, sample


import numpy as np
import pytest

from mmgroup import MM, MMSpace, Cocode
from mmgroup.clifford12 import leech2matrix_eval_A_odd_mod15_aux
from mmgroup.mm import mm_aux_index_extern_to_sparse
from mmgroup.mm import mm_aux_index_leech2_to_sparse
from mmgroup.mm import mm_aux_index_sparse_to_leech2
from mmgroup.mm import mm_aux_index_sparse_to_leech
V = MMSpace(15)



START_AXIS =  ('I', 3, 2)
# The following axes have ben computed externally.
# We have(V(START_AXIS) * MM(AXES[key]).axis_type() == key
# for all kes in the dictionar AXES 
AXES = {
'2A' : 'M<1>',
'4A' : 'M<p_59266661*l_2*p_156183054*l_2*p_54390970*l_1*t_2*y_6c7h*x_651h*d_334h*p_149052354>',
'2B' : 'M<p_133581191*l_1*p_235093658*l_2*p_86199266*l_2*t_1*y_625h*x_9a1h*d_105h*p_200637460>',
'8B' : 'M<p_185765195*l_2*p_72853691*l_2*p_35678269*l_2*t_1*y_10h*x_0dd4h*d_3fh*p_148137261*l_2*p_192868959*l_1*p_17061690*l_1*t_2*y_263h*x_243h*d_206h*p_117493703>',
'6C' : 'M<p_101776500*l_1*p_118232033*l_1*p_236106253*l_1*t_1*y_677h*x_1fc6h*d_353h*p_241875613*l_2*p_165120046*l_2*p_118390844*l_2*t_2*y_3fah*x_9bah*d_64fh*p_216758428>',
'6A' : 'M<p_101776500*l_1*p_118232033*l_1*p_236106253*l_1*t_1*y_677h*x_1fc6h*d_781h*p_4417572*l_2*p_146531284*l_2*p_184837288*l_2*t_1*y_5fah*x_1d1h*d_261h*p_151165491>',
'4C' : 'M<p_29006520*l_2*p_179363783*l_2*p_50934395*l_1*t_1*y_535h*x_14bh*d_298h*p_34465259*l_1*p_213800487*l_1*p_221999319*l_1*t_1*y_5deh*x_0b98h*d_4c8h*p_52448688>',
'4B' : 'M<p_49783959*l_1*p_76678281*l_1*p_187634806*l_1*t_2*y_33ah*x_0e23h*d_258h*p_231417072*l_1*p_145813855*l_1*p_93854609*l_1*t_1*y_2f4h*x_0c20h*d_0deh*p_85411418>',
'6F' : 'M<p_128890944*l_1*p_88397876*l_2*p_29625706*l_1*t_2*y_3e3h*x_0d44h*d_241h*p_72832303*l_2*p_96757696*l_1*p_33880585*l_2*t_1*y_7b8h*x_10eh*d_308h*p_82999734*l_2*p_220813129*l_2*p_161284986*l_2*t_2*y_2abh*x_0d57h*d_1c8h*p_236719595>',
'10B' : 'M<p_79388525*l_1*p_212988967*l_2*p_125281413*l_2*t_2*y_2ebh*x_1a53h*d_6eah*p_243417003*l_2*p_99487229*l_1*p_149175722*l_2*t_2*y_133h*x_723h*d_49fh*p_31110300*l_2*p_225541891*l_2*p_202402963*l_1*t_1*y_2bch*x_75dh*d_4fah*p_65545441>',
'12C' : 'M<p_189176588*l_1*p_49984959*l_1*p_21967454*l_2*t_2*y_7f9h*x_0f72h*d_0b8h*p_83065432*l_1*p_170049576*l_2*p_93174728*l_1*t_2*y_59dh*x_145dh*d_560h*p_100462609*l_1*p_54653496*l_1*p_61568641*l_1*t_2*y_331h*x_0d77h*d_675h*p_219050066>',
'10A' : 'M<p_108457250*l_1*p_152450014*l_1*p_69488708*l_1*t_2*y_74h*x_0c71h*d_368h*p_75469576*l_2*p_21484312*l_2*p_145746011*l_1*t_1*y_242h*x_1e14h*d_671h*p_243860126*l_2*p_147000313*l_1*p_226212428*l_1*t_1*y_439h*x_81ah*d_505h*p_63220017>',
}



BABY_AXES = {
 ('2A', 4):'M<1>',
 ('2A', 0):'M<d_8c9h*p_41220162*l_1*d_720h*p_51535432*l_2*d_9ch*p_122100670*l_1*t_1*y_6eh*x_1d02h*d_31bh>',
 ('2B', 0):'M<d_8c9h*p_41220162*l_1*d_720h*p_51535432*l_2*d_9ch*p_122100670*l_1*t_1*y_6eh*x_1d02h*d_4c0h*p_56421117*l_1*d_794h*p_171856476*l_2*d_616h*p_46300045*l_1*t_2*y_451h*x_19f7h*d_2c5h>',
 ('2B', 8):'M<d_0af8h*p_124682013*l_1*d_373h*p_116663397*l_1*d_732h*p_63395750*l_1*t_2*y_14eh*x_479h*d_0ffch*p_40184923*l_1*d_9c2h*p_60751471*l_2*d_0d5fh*p_43409721*l_1*t_2*y_0eh*x_1015h*d_420h>',
 ('4A', 2):'M<d_8c9h*p_41220162*l_1*d_720h*p_51535432*l_2*d_9ch*p_122100670*l_1*t_1*y_6eh*x_1d02h*d_0e5ch*p_167063353*l_1*d_0ach*p_58851446*l_1*d_6d6h*p_43685151*l_1*t_1*y_5ebh*x_8f8h*d_2deh>',
 ('4B', 1):'M<d_0f04h*p_158228294*l_2*d_9b7h*p_171792597*l_1*d_646h*p_161113830*l_2*t_2*y_4bfh*x_599h*d_476h*p_241545712*l_1*d_18dh*p_175351559*l_2*d_0abh*p_63966998*l_2*t_1*y_104h*x_859h*d_8e1h*p_240348640*l_1*d_206h*p_50612956*l_2*d_0fb3h*p_59975879*l_1*t_2*y_5e6h*x_1052h*d_6dch>',
 ('4C', 2):'M<d_0efdh*p_57294085*l_1*d_82eh*p_152525675*l_2*d_45ch*p_241293137*l_1*t_2*y_449h*x_19e1h*d_436h*p_175340820*l_1*d_155h*p_43760074*l_1*d_6d2h*p_58939823*l_2*t_1*y_576h*x_0c90h*d_382h*p_161111661*l_1*d_1d2h*p_43525311*l_2*d_0a9fh*p_166924760*l_1*t_2*y_10eh*x_18edh*d_6a0h>',
 ('6A', 7):'M<d_0a56h*p_54021241*l_2*d_0d49h*p_122115405*l_2*d_0ebfh*p_54734767*l_2*t_2*y_129h*x_1525h*d_167h*p_49559510*l_2*d_0e5ch*p_159263886*l_1*d_0c95h*p_39263734*l_2*t_2*y_4d2h*x_1df6h*d_59dh*p_161269062*l_1*d_0a73h*p_41274193*l_2*d_442h*p_120011729*l_2*t_2*y_4ddh*x_1467h*d_601h>',
 ('6C', 5):'M<d_184h*p_44294915*l_2*d_71bh*p_128563589*l_1*d_35ch*p_59173956*l_2*t_1*y_169h*x_1de1h*d_296h*p_140266274*l_1*d_686h*p_155930033*l_2*d_0cb2h*p_158249057*l_2*t_1*y_19bh*x_1824h*d_0b2ah*p_52048223*l_2*d_1ebh*p_124185474*l_1*d_514h*p_130867928*l_1*t_2*y_0e5h*x_0ddch*d_7f6h>',
 ('10A', 9):'M<d_0a56h*p_54021241*l_2*d_0d49h*p_122115405*l_2*d_0ebfh*p_54734767*l_2*t_2*y_129h*x_1525h*d_38dh*p_126463070*l_2*d_5afh*p_241453608*l_1*d_0a82h*p_49609083*l_2*t_1*y_1d5h*x_1891h*d_105h*p_59765949*l_2*d_538h*p_164863550*l_2*d_374h*p_134807362*l_2*t_2*y_1d1h*x_0c72h*d_892h*p_241585927*l_2*d_14fh*p_58812896*l_1*d_6a2h*p_53712370*l_2*t_1*y_48eh*x_978h*d_3e9h>',
}


@pytest.mark.involution
def test_axes():
    for key in AXES:
        v = V(START_AXIS) * MM(AXES[key])
        assert v.axis_type() == key
        for i in range(20):
             w = v * MM.rand_G_x0()
             assert w.axis_type() == key
             if (i > 5): continue
             e = randint(1,2)
             we =  w.axis_type(e) 
             assert isinstance(we, str)

