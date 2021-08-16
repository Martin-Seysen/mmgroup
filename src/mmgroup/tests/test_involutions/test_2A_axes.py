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


@pytest.mark.involutions
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

