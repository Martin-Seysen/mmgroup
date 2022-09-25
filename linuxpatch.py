import sys
import os
import subprocess

PATH = 'src/mmgroup'
shared = [x for x in os.listdir(PATH) if x.endswith('.so')] 

def patch_linux():
    if  os.name != 'posix':
        print("This function works for posix systems only")
        return
    ok = True
    for filename in shared:
        path = os.path.join(PATH, filename)
        args = ['patchelf', '--set-path', '$ORIGIN', path]
        try:
            subprocess.check_call(args)
            print(path)
        except:
            print('Executing', ' '.join(args))
            print('failed')
            ok = False
            break
    if ok:
         print("Shared libraries patched successfully")

if __name__ == "__main__":
    patch_linux()    
