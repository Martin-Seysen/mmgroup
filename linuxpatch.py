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



def patch_linux(build_ext_cmd = None, verbose = 1):
    DIR = 'src/mmgroup'
    ABSDIR = os.path.abspath(DIR)
    print([x for x in os.listdir(ABSDIR)])
    shared = [x for x in os.listdir(ABSDIR) if x.endswith('.so')] 
    if  os.name != 'posix':
        print("Warning: Function patch_linux works for posix systems only")
        return
    OUTDIR = None
    if build_ext_cmd:
        build_lib = build_ext_cmd.build_lib
        print("*** build_lib =", build_lib)
        OUTDIR = os.path.abspath(re.sub('src', build_lib, DIR))
    ok = True
    print("Patching shared libraries in Linux...")
    if len(shared) == 0:
        print("Warning: No shared libraries found in directory %s" % ABSDIR)
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
        if OUTDIR:
            dest = os.path.join(OUTDIR, filename)
            if verbose:
                print("Copying %s to %s" % (path, dest))
            copyfile(path, dest)
               
    if ok:
         print("Shared libraries patched successfully")
    else:
         print(MSG)

if __name__ == "__main__":
    patch_linux()

