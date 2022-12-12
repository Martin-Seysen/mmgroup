r"""This little program checks the spider relation in the BiMonster.

"""

try:
    import mmgroup
except (ImportError, ModuleNotFoundError):
    # get the mmgroup package from its inplace location it not installed
    import os
    import sys
    sys.path.append(os.path.join('..', '..', 'src'))
    import mmgroup



from mmgroup.bimm import BiMM, P3_BiMM



def test_spider_relation():
    print("Testing the spider relation in Y_555")
    print("""The spider relation in the BiMonster states that the product
   a * b1 * c1 * a * b2 * c2 * a * b3 * c3
in the BiMonster has order 10.
Here the names of the generators of the BiMonster are as in the Atlas."""
)
    # Enter the spider product into the BiMonster as a string
    spider = P3_BiMM('a,b1,c1,a,b2,c2,a,b3,c3')
    # Check that is has order 10 in the BiMonster
    assert spider.order() == 10, spider.order()
    # Enter the spider product into the BiMonster as a list of generators
    other_spider = ['a', 'b1', 'c1', 'a', 'b2', 'c2', 'a', 'b3', 'c3']
    # Check that its 10th power is equal to the neutral element
    assert P3_BiMM(other_spider * 10) == BiMM(1)
    print("The spider relation has been checked sucessfully")
    




if __name__ == "__main__":
    test_spider_relation()


  