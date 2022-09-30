import sys
import os
import subprocess

PATH = 'src/mmgroup'
shared = [x for x in os.listdir(PATH) if x.endswith('.so')] 


MSG = """
Recommendation:

Install the 'patchelf' utility with the shell command:

sudo apt-get install patchelf

and rebuild the package!

"""



def patch_linux(verbose = 0):
    if  os.name != 'posix':
        print("This function works for posix systems only")
        return
    ok = True
    print("Patching shared libraries in Linux...")
    for filename in shared:
        path = os.path.abspath(os.path.join(PATH, filename))
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
            break
    if ok:
         print("Shared libraries patched successfully")
    else:
         print(MSG)

if __name__ == "__main__":
    patch_linux()

