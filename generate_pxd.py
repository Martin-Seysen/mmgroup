"""Executable program for class generate_code.CodeGenerator"""

import sys




if __name__ == "__main__":
    sys.path.insert(0, "src")
    from mmgroup.generate_c import generate_pxd_parser
    from mmgroup.generate_c import pxdGenerator

    parser = generate_pxd_parser()
    cmdline_args = parser.parse_args(sys.argv[1:])      

    cg = pxdGenerator(cmdline_args)
    cg.activate_py_path()
    cg.import_tables()
    #cg.display_args()
    cg.generate()
    #cg.deactivate_py_path()
    


