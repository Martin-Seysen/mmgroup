from generate_functions import UserDirective
from make_c_tables import TableGenerator 


class TestGen:
    def add_to(x, expr):
         return "{x} = {x} + {expr};\n".format(x = x, expr = str(expr))  

    functions = {
        "AddTo": UserDirective(add_to, "si"),
    }

    tables = {
        "y" : [1,2,3],
        "TABLE_SIZE": 32,
        "mult": 17,
        "pair_list": list(zip(range(3), [7,4,9])),
    }

def generate():
    T = TestGen
    tg = TableGenerator(T.tables, T.functions, verbose = 1)
    tg.generate("testcode.ske", "test.c", "test.h")
    print("Table in test code:")
    print( tg.names )


if __name__ == "__main__":
    generate()

    



