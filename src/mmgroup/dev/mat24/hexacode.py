from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import types

from random import randint

from operator import __xor__
from functools import reduce

def add4(a,b):
    """scalar or vector addition in GF(4)"""
    if type(b) == types.ListType:
         if type(a) == types.ListType:
              assert len(a) == len(b)
              return [a[i] ^ x for i,x in enumerate(b)]
    return a ^ b 


def mul4(a,b):
    """scalar or componentwise vector multiplication in GF(4)"""
    if type(b) == types.ListType:
         if type(a) == types.IntType:
              return [mul4(a,x) for x in b]
         if type(a) == types.ListType:
              assert len(a) == len(b)
              return [mul4(a[i],b[i]) for i in range(len(a))]
    if a * b: return 1 + (a+b-2) % 3
    return 0 


def conj4(a):
    """conjugates a number or vector in GF(4)"""
    if type(a) == types.ListType:
        return [conj4(x) for x in a]
    return [0,1,3,2][a] 
    


def tr4(a):
    """returns trace of a number or componentwise of a vector in GF(4)"""
    if type(a) == types.ListType:
        return [tr4(x) for x in a]
    return a >> 1


def parity(a):
    """returns parity, i.e. sum of entries, of a vector in GF(2) or GF(4)"""
    return reduce(__xor__, a, 0)




def exp4(i):
    """returns \eps**i in GF(4), coded as 2,3,1 for i=0,1,2"""
    return 1 + i%3


def hermite4(a,b):
    """returns tr(a*b**2) for scalars or componentwise for vectors in GF(4)"""
    if type(b) == types.ListType:
         if type(a) == types.IntType:
              return [ hermite4(a,x) for x in b ]
         if type(a) == types.ListType:
              assert len(a) == len(b)
              return [ hermite4(a[i],x) for i,x in enumerate(b) ]
    if a*b: return int(a != b)
    return 0




def quadrilinear_hexacode(a,b,c,d):
   """reuturns scalar product (  tr(a*b**2),  tr(c*d**2) ) """
   h1 = hermite4(a,b)
   h2 = hermite4(c,d)
   return sum([ h1[i] * h2[i] for i in range(6)]) & 1
   

BASIS = [
  [1,0,0,1,3,2],
  [0,1,0,1,2,3],
  [0,0,1,1,1,1]
]


def hexacode_basis(i,j=1):
   """returns j times the i-th hexacode basis vector"""
   return mul4(j, BASIS[i])


def hexacode_random():
   """returns a random hexacode vector"""
   a = [0]*6
   for i in range(3):
      a = add4(a,hexacode_basis(i,randint(0,3)))
   return a


def iter_hexacode():
    """iterator over all 64 hexacode vectors"""
    def mul_basis(i):
        return [ hexacode_basis(i,j) for j in range(4) ]
    for x in mul_basis(0):
        for y in mul_basis(1):
            for z in mul_basis(2):
                 yield(add4(add4(x,y),z))

def iter_hexacode_basis():
    """iterator over nonzero multiples of hexacode basis vectors"""
    for i in range(3):
        for j in range(1,4):
            yield hexacode_basis(i,j)

def weight(a):
    """returns the weight of a vector. i.e. the No of nonzero entries"""
    return sum(map(bool,a))


def cartesian(*lists):
    """returns the Cartesian product of the argument, each argument a list"""
    if len(lists) == 0:
        yield []
        return
    hd, tl = lists[0], lists[1:]
    for h in hd:
        for l in cartesian(*tl):
            yield [h] + l





def cocode_product(a,b):
    c = conj4(mul4(a,b))
    for i in range(3):
        c = add4(c, hexacode_basis(i,c[i]) )
    for i in range(3): assert c[i] == 0
    return c 


def reduce_to_cocode(a):
    c = a[:]
    for i in range(3):
        c = add4(c, hexacode_basis(i,c[i]) )
    for i in range(3): assert c[i] == 0
    return c 

def is_hexacode_word(a):
    return reduce_to_cocode(a) == [0] * 6
  
LIGHT_COCODE_WORDS = None

def light_cocode_words():
    """returns dictionary mapping cocode words to lists of light cocode words

    Here the keys in the dictionary are the possible return values of
    function reduce_to_cocode().

    The light cocode words are the words of weight <= 2
    """
    ## TODO: cleanup and document this !!!!!!!!!!!!!!!!
    global LIGHT_COCODE_WORDS
    if not LIGHT_COCODE_WORDS is None:
        return LIGHT_COCODE_WORDS
    LIGHT_COCODE_WORDS = {}
    def add_word(cocode_word, light_word):
        cw = tuple(cocode_word)
        if not cw in LIGHT_COCODE_WORDS.keys():
             LIGHT_COCODE_WORDS[cw] = [ light_word[:] ]
        else:
             LIGHT_COCODE_WORDS[cw] = ( LIGHT_COCODE_WORDS[cw] + [light_word[:]] )
    add_word([0]*6, [0]*6) 
    for i in range(6):
         for wi in range(1,4):
             lw = [0]*i + [wi] +  [0]*(5-i)
             cw = reduce_to_cocode(lw)
             add_word(cw, lw)
             for j in range(i+1,6):
                 for wj in range(1,4):
                     lw2 = lw[:]
                     lw2[j] = wj
                     cw = reduce_to_cocode(lw2)
                     add_word(cw, lw2)
    return LIGHT_COCODE_WORDS 


def decode_cocode(a):
    """returns list of equivalent words to a of length <= 2 modulo the hexacode"""
    cw = tuple(reduce_to_cocode(a))  
    return  light_cocode_words()[cw]  

def codeword_type(a):
    """returns a tuple describing the type of a vector a in GF(4)^6.

    The 1st entry is the weight of the vector a. Other entries are present
    if a is no hexacode word. Then the other entries are the (sorted) tuple 
    of the weights of all hexacode words that have minimal  distance to a.

    """
    ## Baustelle: cleanup and document this !!!!!!!!!!!!!!!!
    wa = weight(a)
    wd = decode_cocode(a)
    if weight(wd[0]) == 0: return (wa,)
    cwd = [ add4(a,b) for b in wd]
    wd1 = sorted(map(weight, cwd))
    return tuple( [wa]  +  wd1 )

            


subword4_list = None

def subwords_4():
    """return pattern list of all words of weight 4"""
    global subword4_list
    if subword4_list is None:
         subword4_list = []
         for i in range(6):
             for j in range(i+1,6):
                 subword4_list.append( 
                     [3]*(i) + [0] + [3]*(j-i-1) + [0] + [3]*(6-j-1) )
    return subword4_list

def no_hexacode_subwords(a):
    """Returns True if no subword of a is a hexacode word

    subwords are obtained by replacing entries of a by zero.
    """
    if is_hexacode_word(a):
        return False
    for y in subwords_4():
        b = [a[i] & y[i] for i in range(6)]
        if is_hexacode_word(b):
            return False
    return True
 

def iter_content(length=6, content = [1,2,3]):
    """Iterates over all lists of given length with given content"""
    if length == 0:
        yield []
    else:
        for l in iter_content(length-1, content):
            for j in content:
                 yield l + [j]


    
  


def trace3(a,b,c):
    """return trace(a*b*c), a,b,c \in GF(4)"""
    return mul4(a,mul4(b,c)) >> 1

def trilinear(a,b,c):
    
    y = 0
    for i in range(6):
        y ^= trace3(a[i],b[i],c[i])
    return y 





###############################################################################




def test_hermite4():
    for a in range(4):
        for b in range(4):
            assert hermite4(a,b) == tr4(mul4(a, conj4(b))) 


def test_quadrilinear_nontrivial():
    print( "\nTest quadrilinear form on QA random hexacode vectors .." )
    for i in range(100):
        a = hexacode_random()
        b = hexacode_random()
        c = hexacode_random()
        d = hexacode_random()
        x = quadrilinear_hexacode(a,b,c,d) 
        y = quadrilinear_hexacode(a,c,b,d) 
        z = quadrilinear_hexacode(a,d,c,b)
        if x == y == z == 0: continue  
        print( a, b, c, d )
        print( x,y,z )
        print( "quadrilinear form QA is not trivial\n" )
        return



def test_quadrilinear(a,b,c,d):
    A = quadrilinear_hexacode(a,b,c,d) 
    B = quadrilinear_hexacode(a,c,b,d) 
    C = quadrilinear_hexacode(a,d,c,b) 
    assert A ^ B ^ C == 0 


def test_quadrilinear_on_basis():
    for a in iter_hexacode_basis():
        for  b in iter_hexacode_basis():
             for  c in iter_hexacode_basis():
                 for  d in iter_hexacode_basis():
                     test_quadrilinear(a,b,c,d)
    print( "QA(a,b,c,d) + QA(a,c,b,d) + QA(a,d,c,b) = 0""" )
    for a in iter_hexacode():
        for b in iter_hexacode():
            assert parity( hermite4(a,b) ) == 0   
    print( "QA(a,b,a,b) = 0" )




def test_trilinear_form():
    print( "put A(a,b,c) = sum_{i} tr(a[i] b[i] c[i])" )
    for a in iter_hexacode_basis(): 
        for b in iter_hexacode_basis():
            for c in iter_hexacode_basis():
                if  trilinear(a,b,c):
                    print( "A(",a,',',b,',',c,") = 1" )
                    print( "trilinear form A : hexacode -> GF(2) is not trivial" )
                    return

def test_self_dual():
    for a in iter_hexacode(): 
        for b in iter_hexacode():
            assert sum(hermite4(a,b)) & 1 == 0
    print( "Hexacode is self-dual, i.e. hermite4(a,b) = 0 for codewords a,b" )


def hermite4_cocode_nonzero():
   print( "Example where hermite4(a,b) = (..,tr(a_i*b_i^2),..) is zero" )
   print( "and cocode product A(a,b,.) is not zero" )
   null = [0]*6
   for a in iter_hexacode(): 
      for b in iter_hexacode():
          c = cocode_product(a,b)
          if  c != null:
              try:  
                  assert hermite4(a,b) != null, (a,b,c, hermite4(a,b) )
              except: 
                  print( "hermite4(",a,',',b,") = 0" )
                  print( "cocode product is",  conj4(mul4(a,b)) )  
                  return




def print_no_hexacode_subwords():
    n = 0
    for a in iter_content(6,  content = [1,2] ):
        if  no_hexacode_subwords(a) and sum(a) == 8:
            n += 1
    print ("\nThere are", n, "words containing four times 1 and twice 2")

    L2 = []
    print( " such that no subword is a hexacode word" )
    for x0 in iter_content(4,  content = [2,3] ):
       x = [1,1] + x0
       if sum(x) == 12 and  no_hexacode_subwords(x):
           L2.append(x)
    print( "For the following words no subword is a hexacode word:" )
    print( L2 )




def hexacode_types():
    # We go thru the Cartesian product [0,1,2,3] x ... x [0,1,2,3] with 6 factors"""
    all_words = list(cartesian( *([[0,1,2,3]]*6)))
    print("""\nDetermine types of all words of GF(4)^6 with respect to the hexacode:

The 1st entry is the weight of the vector v. Other entries are present
if v is no hexacode word. Then the other entries are the (sorted) tuple 
of the weights of all hexacode words that have minimal distance to v.""")
    types = {}
    samples = {}
    for a in all_words:
        newtype = codeword_type(a)
        if not newtype in types:
             types[newtype] =0
             samples[newtype] = (a, decode_cocode(a))
        types[newtype] += 1
    n = 0
    for k in sorted(types.keys()):
         if  k[0] != n: print( '' )
         print( str(k) + ':' + str(types[k]) +',', end=' ' )
         n =k[0] 
    return types
 
    




#######################################################################



def selftest():
    s ="""Let H be the Hexacode over GF(4)**6, QA: H**4 -> GF(2) the quadrilinear
form given by QA(a,b,c,d) = (  tr(a*b**2),  tr(c*d**2) ) with (.,.) the 
scalar product on GF(2)**6 and tr(a*b**2) computed componentwise on the 
hexcode.  Note that  tr(a*b**2) is in GF(2) for a,b in GF(4),
"""

    print( s )

    print("Hexacode basis vectors and their multiples, with '2' a generator of GF(4):")
    for i in range(3):
       for j in range(1,4):
         print ( hexacode_basis(i,j), ',' , end=' ' )
       print ("")


    test_quadrilinear_nontrivial()

    test_quadrilinear_on_basis()

    print ("")

    test_trilinear_form()
    print ("")
                   
    test_hermite4()

    hermite4_cocode_nonzero()

    test_self_dual()


    print_no_hexacode_subwords()


    hexacode_types()


if __name__ == "__main__":
    selftest()

