from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import zipfile, time
import os
import subprocess




try:
    import zlib
    compression = zipfile.ZIP_DEFLATED
except:
    compression = zipfile.ZIP_STORED
    raise

modes = { zipfile.ZIP_DEFLATED: 'deflated',
          zipfile.ZIP_STORED:   'stored',
          }

excluded_dirs = [ 'backup', 'eggs', 'EGG-INFO', '__pycache__', 
    'build', 'c_files', 'c_doc',
    '_doc', '.pytest_cache', 'immaturre_tests', ]

included_extensions = [".bat"]
excluded_extensions = [
    ".pyd", ".pyc", ".o", ".bak", ".s", ".log", ".def", ".dll", 
]

def zipdir(path, zip):
    for root, dirs, files in os.walk(path):  
        root_tail = os.path.split(root)[1]
        print( os.path.split(root), dirs)
        if root_tail in excluded_dirs:
            continue     
        for file in files:
            _ , ext = os.path.splitext(file) 
            if ext in included_extensions or not root in excluded_dirs:
                 if ext not in excluded_extensions:
                     pathname = os.path.join(root, file)                     
                     print(pathname)
                     zip.write(pathname, compress_type=compression)



tm = time.strftime("%Y_%m_%d_%H_%M")

filename = r'Monster_%s.zip' % tm

zip_path = os.path.join('backup', filename)

backup_dir = r"E:\AlterPC\projects\MonsterNew\backup"


if __name__ == '__main__':
    zipf = zipfile.ZipFile(zip_path, 'x')
    zipdir('.', zipf)
    zipf.close()
    try:
        print( subprocess.check_output(["cmd", "/c", "copy", zip_path, backup_dir]) )
        zf = zipfile.ZipFile(os.path.join(backup_dir, os.path.join(backup_dir,filename)))
        print ("Backup to external drive successful!")
    except:
        print ("Backup to external drive did not work!")
        raise

