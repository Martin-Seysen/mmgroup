r"""This little program checks the spider relation in the Bimonster.

"""

try:
    import mmgroup
except (ImportError, ModuleNotFoundError):
    # get the mmgroup package from its inplace location if not installed
    import os
    import sys
    sys.path.append(os.path.join('..', '..', 'src'))
    import mmgroup



from mmgroup.bimm import BiMM, P3_BiMM



def test_spider_relation():
    print("Testing the spider relation in Y_555")
    print("""The spider relation in the Bimonster states that the product
   a * b_1 * c_1 * a * b_2 * c_2 * a * b_3 * c_3
in the Bimonster has order 10.
Here the names of the generators of the Bimonster are as in the Atlas."""
)

    # Here is one test of the spider relation
    # List of generators of IncP3 corresponding to the spider relation 
    spider = ['a', 'b1', 'c1', 'a', 'b2', 'c2', 'a', 'b3', 'c3'] * 10
    # Let 'Spider' be the image of the spider relation in the Bimonster
    Spider = P3_BiMM(spider)
    # Check that this is the neutral element of the Bimonster
    assert Spider == BiMM(1)

    # Here is another test of the spider relation
    # Put Spider1 = a * b_1 * c_1 * a * b_2 * c_2 * a * b_3 * c_3.
    # Here we enter that product into the Bimonster as a string
    Spider1 = P3_BiMM('a, b1, c1, a, b2, c2, a, b3, c3')
    # Check that Spider1 has order 10 in the Bimonster
    assert Spider1.order() == 10

    print("The spider relation has been checked sucessfully")
    




if __name__ == "__main__":
    test_spider_relation()


  