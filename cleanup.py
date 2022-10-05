import os
import subprocess
from optparse import OptionParser

parser = OptionParser()


ROOT_DIR = os.path.dirname(os.path.realpath(__file__))


PROGRAM_FILES = [".c", ".h", ".pxd", ".pxi", ".pyx"]

EXTENSIONS = [".pyd", ".dll", ".lib", ".so"]


EXTENSION_STR = ", ".join(EXTENSIONS) if len(EXTENSIONS) else "?"



def del_pyc(verbose = False):
    if verbose:
        print("\nDeleting python intermedate files (.pyc)") 
    for root, dirs, files in os.walk(ROOT_DIR, topdown=False):
        for name in files:
            if os.path.splitext(name)[1] == ".pyc":
                path = os.path.join(root, name)
                try:
                    if verbose:
                       print("Delete %s" % path)
                    os.remove(path)
                except:
                    if verbose:
                       print("failed")


def del_c(verbose = False):
    C_DIRS = [
        ["src", "mmgroup", "dev", "c_files"],
        ["src", "mmgroup", "dev", "pxd_files"],
    ]
    if verbose:
        print("\nDeleting automatically generated .c and .pxd files") 
    for dir_list in C_DIRS:
        try:
            dir = os.path.join(ROOT_DIR, *dir_list)
        except:
            continue
        files = os.listdir(dir)
        for name in files:
            if name == "readme.txt":
                continue
            path = os.path.join(dir, name)
            try:
                if verbose:
                   print("Delete %s" % path)
                os.remove(path)
            except:
                if verbose:
                   print("failed")


DATA_FILE_DICT = {
    ("src", "mmgroup", "dev", "mm_reduce" ): [
        "order_vector_data", "v1_mod3_data", 
    ],
    ("src", "mmgroup", "tests", "test_axes"): [
        "sample_axes", "baby_sample_axes",
    ],
    ("src", "mmgroup", "tests", "test_involutions"): [
        "involution_samples",
    ]
}



def del_data(verbose = False):
    if verbose:
        print("\nDeleting automatically generated python data files") 
    for path_info, file_list in DATA_FILE_DICT.items():
        dir = os.path.join(ROOT_DIR, *path_info)
        for filename in file_list:
            path = os.path.join(dir, filename + ".py")
            try:
                if verbose:
                   print("Delete %s" % path)
                os.remove(path)
            except:
                if verbose:
                   print("failed")

def del_ext(verbose = False):
    if verbose:
        print("\nDeleting automatically generated python extensions") 
    EXT_DIRS = [
        ["src", "mmgroup", "dev", "c_files"],
        ["src", "mmgroup"],
    ]
    if len(EXTENSIONS):
        for ext_dir in EXT_DIRS:
            dir = os.path.join(ROOT_DIR, *ext_dir)
            files = os.listdir(dir)
            for name in files:
                if os.path.splitext(name)[1] in EXTENSIONS:
                    path = os.path.join(dir, name)
                    try:
                        if verbose:
                           print("Delete %s" % path)
                        os.remove(path)
                    except:
                        if verbose:
                            print("failed")
    else:
        s = "Dont't know how do delete python extensions in '%s' system"
        print(s % os.name)
      

def git_checkout_data(verbose = False):
    if verbose:
        print("\nCheck out automatically generated data files with git") 
    for path_info, file_list in DATA_FILE_DICT.items():
        dir = os.path.join(*path_info)
        for filename in file_list:
            path = os.path.join(dir, filename + ".py")
            try:
                if verbose:
                   print("Checking out %s" % path)
                args = ["git", "checkout", path]
                print(" ".join(args))
                subprocess.check_call(args)
            except:
                if verbose:
                   print("failed")




def main():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option("-p",  dest="del_pyc", action="store_true",
        help="Delete intermediate python files (.pyc)")
    parser.add_option("-c",  dest="del_c", action="store_true",
        help="Delete automatically generated program files (%s)" %
        ", ".join(PROGRAM_FILES))
    parser.add_option("-d",  dest="del_data", action="store_true",
        help="Delete automatically generated data files (.py)")
    parser.add_option("-x",  dest="del_ext", action="store_true",
        help="Delete automatically generated extensions (%s)" %
        EXTENSION_STR)
    parser.add_option("-a",  dest="del_all", action="store_true",
        help="Delete all automatically generated files" )
    parser.add_option("-g",  dest="git_checkout", action="store_true",
        help="Checkout all automatically generated files with git" )
    parser.add_option("-v",  dest="verbose", action="store_true",
        help="Verbose operation" )
    
    options, args = parser.parse_args()
    verbose = options.verbose

    if options.del_pyc or options.del_all:
        del_pyc(verbose)
    if options.del_c or options.del_all:
        del_c(verbose)
    if options.del_data or options.del_all:
        del_data(verbose)
    if options.del_ext or options.del_all:
        del_ext(verbose)
    if options.git_checkout:
        git_checkout_data(verbose)


if __name__ == "__main__":
    main()


