import sys
import os
import subprocess
import  re
from shutil import copyfile



MSG = """
Recommendation:

Install the 'patchelf' utility with the shell command:

sudo apt-get install patchelf

and rebuild the package!

"""


def has_extension(filename, extension_list):
    return any(filename.endswith(ext) for ext in extension_list)



def patch_and_copy_shared(build_ext_cmd = None, verbose = 1):
    DIR = 'src/mmgroup'
    ABSDIR = os.path.abspath(DIR)
    if os.name in ["nt"]:
        extensions =  [".pyd", ".dll"]     
    elif os.name in ["posix"]:
        extensions = [".so"]
    else:  
        W = "don't know how do process shared libraries in a %s system"
        print("Warning:", W % os.name)     
        extensions =  []  

    files = os.listdir(ABSDIR)
    #print(files)
    shared = [x for x in files if has_extension(x, extensions)] 
    OUTDIR = None
    if build_ext_cmd and not build_ext_cmd.inplace:
        build_lib = build_ext_cmd.build_lib
        print("*** build_lib in setup.py =", build_lib)
        OUTDIR = os.path.abspath(DIR.replace('src', build_lib))
    if len(shared) == 0:
        print("Warning: No shared libraries found in directory %s" % ABSDIR)

    if os.name in ["posix"]:
        ok = True
        print("Patching shared libraries in Linux...")
        for filename in shared:
            path = os.path.abspath(os.path.join(ABSDIR, filename))
            if verbose:
                print("patching " + path)
            args = ['patchelf', '--set-rpath', '$ORIGIN', path]
            try:
                subprocess.check_call(args)
                #print(path)
            except:
                print('Executing', ' '.join(args))
                print('failed')
                ok = False
        if ok:
            print("Shared libraries patched successfully")
        else:
            print(MSG)

    if OUTDIR:
        for filename in shared:
            path = os.path.abspath(os.path.join(ABSDIR, filename))
            dest = os.path.join(OUTDIR, filename)
            if verbose:
                print("Copying %s to %s" % (path, dest))
            copyfile(path, dest)
               

if __name__ == "__main__":
    patch_and_copy_shared()

