"""Executable program for class generate_code.CodeGenerator"""

import sys




if __name__ == "__main__":
    sys.path.insert(0, "src")
    from mmgroup.generate_c import generate_code_parser
    from mmgroup.generate_c import CodeGenerator

    parser = generate_code_parser()
    cmdline_args = parser.parse_args(sys.argv[1:])      
    #print("Arguments:", sys.argv)

    cg = CodeGenerator(cmdline_args)

    cg.activate_py_path()
    cg.import_tables()
    cg.display_args()
    cg.generate()
    #cg.deactivate_py_path()
    


