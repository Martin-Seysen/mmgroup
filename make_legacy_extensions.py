import sys
import argparse




class MyArgumentParser(argparse.ArgumentParser):
    def convert_arg_line_to_args(self, arg_line):
        return arg_line.split()




def generate_legacy_code_parser():
    description = ("Generate python substitutes for legacy entensions"
    )

    parser = MyArgumentParser(
        description=description
    )

    parser.add_argument('--out-dir', 
        metavar='OUT_DIR',
        action = 'store',
        help="Generate python substitutes for legacy entensions"
             "in directory OUT_DIR"
    )
    return parser



if __name__ == "__main__":
    sys.path.insert(0, "src")
    from mmgroup.dev.mm_op.dispatch_p import make_all_legacy_scripts
    parser = generate_legacy_code_parser()
    cmdline_args = parser.parse_args(sys.argv[1:])      
    make_all_legacy_scripts(cmdline_args.out_dir)

