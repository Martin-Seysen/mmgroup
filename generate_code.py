"""Executable program for class generate_code.CodeGenerator"""

import sys




if __name__ == "__main__":
    sys.path.insert(0, "src")
    from mmgroup.generate_c import generate_code_parser
    from mmgroup.generate_c import CodeGenerator

    parser = generate_code_parser()
    parser.add_argument('--output-c-names', 
        default = False,
        action = 'store_true',
        help = 'Output names of C files to be generated and exit'
    )
    cmdline_args = parser.parse_args(sys.argv[1:])      
    #print("Arguments:", sys.argv)

    cg = CodeGenerator(cmdline_args)
    if cmdline_args.output_c_names:
        c_files = cg.c_files()
        for file in c_files:
            print(file)
        sys.exit(0)

    cg.activate_py_path()
    cg.import_tables()
    cg.display_args()
    cg.generate()
    #cg.deactivate_py_path()
    


