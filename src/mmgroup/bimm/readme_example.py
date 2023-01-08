r"""
As stated in the description of the group :math:`Y_{555}` , the spider relation 
in that group is:

.. math:: 

   (a b_1 c_1 a b_2 c_2 a b_3 c_3)^{10} = 1 \, . 

We can quickly check the spider relation as follows:

.. code-block:: python

    # Class BiMM implements an element of the Bimonster.
    # Function P3_BiMM maps a product of generators of 
    # the Coxeter group IncP3 to the Bimonster.
    from mmgroup.bimm import  BiMM, P3_BiMM
    # List of generators of IncP3 corresponding to the spider relation 
    spider = ['a', 'b1', 'c1', 'a', 'b2', 'c2', 'a', 'b3', 'c3'] * 10
    # Let 'Spider' be the image of the spider relation in the Bimonster
    Spider = P3_BiMM(spider)
    # Check that this is equal to >the neutral element of the Bimonster
    assert Spider == BiMM(1)


Alternatively, the spider relation can be verified as follows:

.. code-block:: python

    from mmgroup.bimm import  P3_BiMM
    # Put Spider1 = a * b_1 * c_1 * a * b_2 * c_2 * a * b_3 * c_3.
    # Here we enter that product into the Bimonster as a string
    Spider1 = P3_BiMM('a, b1, c1, a, b2, c2, a, b3, c3')
    # Check that Spider1 has order 10 in the Bimonster
    assert Spider1.order() == 10




"""

