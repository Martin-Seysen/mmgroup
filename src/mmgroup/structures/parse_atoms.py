"""Parser for obtaining atoms of algebraic structures

This module contains classes an functions that support parsing of 
atoms of an algebraic structure from a string. In this module we
abuse the word "algebra" for denoting any algebraic structure.

Here an algebra is represented as subclass of an abstract class 
modelling an abstract algebra. Examples of such abstract classes 
are AbstractGroup modelling a group, and AbstractRepSpace modelling
a free finite-dimensional module over the integers modulo p with
a group operating on that module.

In order to make an abstract class AbstractXxx useful there are 
subclasses of class AbstractXxx, such that an instance of such a
subclass models a specific algebraic structure. Such a subclass
contains functions for performing the operations on the algebraic
structure. For each class AbstractXxx there is a class 
AbstactXxxWord modelling an element of a specific algebra. For 
an instance x of AbstactXxxWord the function x.algebra() returns 
the specific algebra to which x belongs. The main purpose of the 
class AbstactXxxWord is to overload standard python operators with 
the appropriate functions defined in the corresponding algebra. 
Since each class modelling an algebra contains functions to 
generate elements of that algebra, users need not deal with class 
AbstactXxxWord.


Many algebras are built up from atoms, e.g. groups may be built
up from generators and vector spaces may be built up from
multiples of unit vectors. There is a standard method to map
some strings to atoms for a specific algebra.

Here all atoms are tuples of the form

    (tag, data1, data2, ...)

where 'tag' is an alphabetic string of length 1 and 
data1, data2, ... are either nonnegative integers or alphanumeric
strings representing identifiers and not containing the undersore
character '_'.

We support a parser that reckognizes as string of shape

    tag_data1_data2[_ ...]

as an atom. Such a parser can distinguish between atoms, idenifiers
and numbers, where numbers are nonnegative integers.

The grammar for distingushing between atoms, identifiers and numbers 
in Extended Backus-Naur form is as follows:


literal           =  atom | identifer | number ;
atom              =  tag, "_" , index, {"_",  index} ;
identifier        =  (letter, (letter | digit), {alphanum}) ; 
number            =  decimal number | hex number ; 
tag               =  letter ;
index             =  num index | alpha index ;
num index         =  decimal number | hex index ;
hex index         =  hex number | (digit, {hex digit}, ("h" | "H")) ;
decimal number    =  digit, {digit} ;
hex number        =  "0x", hex digit, {hex digit} ;
alpha index       =  letter, {letter | digit} ;
alphanum          =  letter | digit | "_" ;
hex digit         =  digit | "a" ... "f" | "A" ... "F" ;
letter            =  "A" ... "Z" | "a" ... "z" ;
digit             =  "0" ... "9" ;

It turns out that the standard python parser reckognizes both,
atoms and identifiers as python identifers. 

Given a specific algebra, we may evaluate any python expression 
(using e.g. the python eval function) in such a way that all
python identifiers representing atoms are evaluated to the 
corresponding atoms of the specific algebra. Therefore we must
pass a mapping from the identifier strings to their values to
the evaluation function. Then the evaluation function parses
a string as a python expression, it maps all identifiers to
their value an returns the evaluated expression.

For an abstract algebra class AbstractXxx there may be an abstract 
class AtomicXxx modelling an algebra of type AbstractXxx generated
by atoms. We call such an algebra an 'atomic' algebra. As in the
case of abstract algebras, a specific algebra in modelled as an
instance of a subclass of class AbstractXxx. 

Calling a specific algebra with a single argument converts this
argument to an element of the algebra, as usual in python. If
the specific algebra is an atomic algebra, that parameter may be 
a python expression coded as a string. The idenifiers in that 
string, which have the shape of an atom according to the grammar
given above, are interpreted as atoms. This provides a simple
mechanism for generating atoms of an algebra.
 
"""

import re
import ast
import warnings


#######################################################################
# Class AtomDict
#######################################################################


class AtomDict(object):
    """Models a parser for obtaining atoms as described above

    The constructor of this class takes two arguments :
    - a function 'reducer', 
    - a dictionary 'identifiers'.

    For an instance a of this class and a string s, a[s]
    parses the string a according to the grammer given above and
    returns one of the following objects:

    - an integer
    - the result of a function f(tag, data1, data2,...)
    - an object referred by a valid python identifier.

    If s is an integer recognized by the grammer given above
    then a[s] returns that integer.

    If s is a string of shape  tag_data1_data2[_ ...] as described
    by the grammar given above, then a[s] returns 
    reducer(tag, data1, data2, ...). Here 'tag' is a letter. i.e. a
    string of length 1, and the substrings data1, data2, ...  are
    evaluated to strings or to integers  according to the grammar
    given above.

    If s is an identifier contained as a key in the dictionary
    'identifiers' then a[s] returns identifiers[s].

    In any other case the invocation of a[s] raises ValueError.
    """
    match_atom = re.compile(r"[A-Za-z]_\w+", re.ASCII)
    match_identifier = re.compile(r"[A-Za-z][A-Za-z0-9]\w+", re.ASCII)
    match_number = re.compile(r"[0-9]+|0x[0-9A-Fa-f]+")
    match_letter = re.compile(r"[A-Za-z]")

    m_int = re.compile(r"([0-9]+)")
    m_hex1 = re.compile(r"([0-9][0-9A-Fa-f]*)[Hh]")
    m_hex2 = re.compile(r"(0x[0-9A-Fa-f]+)")
    m_random = re.compile(r"([A-Za-z][A-Za-z0-9]*)")

    m_dict_index = {
        m_int: int,
        m_hex1: lambda x : int(x,16),
        m_hex2: eval,
        m_random: str,
    }

    def __init__(self, reducer = tuple, identifiers = {}):
        self.identifiers = identifiers
        self.reducer = reducer

    def __getitem__(self, s):
        if self.match_atom.fullmatch(s):
            indices = map(self._parse_index, s[2:].split("_"))
            word = self.reducer(s[0], *indices)
            return word
        elif self.match_identifier.fullmatch(s):
            try:
                return self.identifiers[s]
            except KeyError:
                raise KeyError("Identifer %s not found in atom literal" %s)
        elif self.match_number.fullmatch(s):
            return eval(s)
        raise ValueError("Syntax error in atom literal")

    def _parse_index(self, s):
        for matobj, translator in self.m_dict_index.items():
            try:
                return translator(matobj.fullmatch(s).groups()[0])
            except:
                pass
        return ValueError("Syntax error in index of atom literal")

       


#######################################################################
# Evaluating a string with function eval_atom_expression()
#######################################################################


# For background on module ast (abstact syntax tree) see e.g.:
# https://greentreesnakes.readthedocs.io/en/latest/tofrom.html







# Set of nodes in module .ast considered as safe for function eval()
# See  https://greentreesnakes.readthedocs.io/en/latest/nodes.html
good_ast_nodes = set([
   ast.Expression,
   ast.Expr,
   ast.Load,
   ast.Constant,

   ast.UnaryOp,
   ast.UAdd, ast.USub, ast.Not, ast.Invert,

   ast.BinOp,
   ast.Add, ast.Sub, ast.Mult, ast.Div, ast.FloorDiv, ast.Mod, 
   ast.Pow, ast.LShift, ast.RShift, ast.BitOr, ast.BitXor, 
   ast.BitAnd,

   ast.Call,
   ast.Subscript, ast.Index 
])


with warnings.catch_warnings():
    warnings.filterwarnings('error')
    # Some deprecated ast classes
    AST_CLASSES = ['Num', 'Str', 'NameConstant']
    for attr in AST_CLASSES:
        try:
            good_ast_nodes.add(getattr(ast, attr))
        except:
            pass



class EvalNodeVisitor(ast.NodeVisitor):
    """Node visitor for Abstract Syntax Trees creeated by module  .ast 

    This vistor checks for all nodes of an AST whether they are secure 
    for processing with function eval(). Method visit() raises TypeError
    if an insecure node id found in the AST.

    Caution: 
    Leading underscores in python identifiers or attributes are illegal! 

    Note that leading underscores in an identifier or an attribute are
    considered dangerous. E.g.  parsing '().__class__.__bases__[0]'
    yields the class <object>, which a user should not manipulate. For  
    a discussion about the security of eval() see:
    https://nedbatchelder.com/blog/201206/eval_really_is_dangerous.html

    Creating tuples in backets or dictionaries in curly brackets is also
    illegal because the syntax gets messy when used in the code generator.
    """
    def visit(self, Node):
        cls = Node.__class__
        #print("visiting", Node, cls)
        if cls in good_ast_nodes:
            pass
        elif cls == ast.Name:
            if Node.id[:1] == "_":
                raise NameError("Illegal leading '_' in name in python expression") 
        elif cls == ast.Attribute:
            if Node.attr[:1] == "_":
                raise NameError("Illegal leading '_' in name in python expression") 
        else:
            print("\nIllegal ast node", type(Node))
            raise TypeError("Illegal construct in python expression") 
        self.generic_visit(Node)


def eval_atom_expression(s, atom_dict):
    """An attempt to create a secure version of function eval()

    See function EvalNodeVisitor() for allowd python contructs.
    """
    try:
        a_tree = ast.parse(s.strip(), mode="eval")
        EvalNodeVisitor().visit(a_tree)
        code = compile(a_tree, filename="atom expression", mode="eval")
        result = eval(code, {'__builtins__': {}},  atom_dict )
        return result
    except:
        print("\nError in evaluating an atom expression!\n")
        raise
        

#######################################################################
# Class TaggedAtom
#######################################################################

def ihex(x, length = None):
    x = int(x)
    sign = "-"[:(x < 0)]
    s = "%xh" % abs(x) if length is None else "%0*xh" % (length, abs(x))
    return sign + s if s[:1].isdigit() else sign + "0" + s

class TaggedAtom(object):
    """Models an atom of a general group with generators and relations"""

    # format_tuple[i] is used for formatting the i-th component.
    # Recommended entries of the tuple are str, hex, and ihex
    format_tuple = (str,)
    operator_char = "*"
      
    def __init__(self, tag, *data):
        """Construct an atom from a tag and a seqence of od data. 

        'tag' should be a single alphabetical character. All
        data should be nonnegative integers or alphanumeric
        strings with the first character not a digit.

        Subclasses may overwerite this method.
        """
        self.data = data
        assert isinstance(tag, str) and len(tag) == 1
        self.tag = tag
        
    def as_tuples(self):
        """Return representation of an atom as a list of tuples

        Subclasses may overwerite this method. The function
        may be implemented as a generator
        """
        yield (self.tag,) + tuple(int(x) for x in self.data)

    def reduce(self, *args):
        """Reduce atom 'self' to a simpler form if possible"""
        return self
    	
    def str_tuple(self, tuple_):
        strings = []
        for i, x in enumerate(tuple_):
            try:
                f = self.format_tuple[i]
            except IndexError:
                f = str
            strings.append(f(x))
        return "_".join(strings)

    def str(self):
        tuples = []
        for tuple_ in self.as_tuples():
            tuples.append(self.str_tuple(tuple_))  
        return self.operator_char.join(tuples)   

    __repr__ = str
    



#######################################################################
# Test program
#######################################################################


if __name__ == "__main__":
    def atom_to_tuple(*args):
        yield args
    def atom_to_part1(*args):
        yield args[1]
    def reducer(*data):
       return data[0]
    atom_dict = {"t": atom_to_tuple, "u" : atom_to_part1}
    a = AtomDict(reducer, atom_dict, {"x17": 17.5})
    print( a["t_1_a"], a["x17"])
    print( eval_atom_expression("(3 + u_1 ) * 100 * x17", a))
        