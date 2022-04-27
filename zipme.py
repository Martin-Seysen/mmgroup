"""Zip the whole project to a zip file in subdirectory 'build'.


The zip process is based on:

python setup.py sdist --format=zip

Non-source files are removed from the zip file.

"""


from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import zipfile, time
import sys, os
import subprocess

from zipfile import ZIP_DEFLATED


tm = time.strftime("%Y_%m_%d_%H_%M")

filename = r'Monster_%s.zip' % tm


source = r"dist\mmgroup-0.0.7.zip"

zip_path = os.path.join('backup', filename)

backup_dir = r"E:\AlterPC\projects\MonsterSup\backup"




def is_to_be_zipped(filename):
    if filename[-4:] in [".pyd", ".dll"]:
        return False
    if filename.find("build/") >= 0:
        return False 
    if filename.find("docs/build/") >= 0:
        return False 
    if filename.find("egg-info") >= 0:
        return False 
    if filename.find("pytest_cache") >= 0:
        return False 
    if filename.find("/src/mmgroup/dev/c_files/") >= 0:
        return False 
    if filename.find("/src/mmgroup/dev/c_doc/") >= 0:
        return False 
    if filename.find("/src/mmgroup/dev/pxd_files/") >= 0:
        return False        
    return True


if __name__ == '__main__':
    if os.path.exists(source):
        os.remove(source)

    subprocess.check_output([sys.executable, "setup.py", "sdist", "--format=zip"]) 

    zin = zipfile.ZipFile (source, 'r')
    zout = zipfile.ZipFile (zip_path, 'w', compression=ZIP_DEFLATED,  compresslevel=9)
    for item in zin.infolist():
        buffer = zin.read(item.filename)
        #print(item.filename, is_to_be_zipped(item.filename))
        if is_to_be_zipped(item.filename):
            zout.writestr(item, buffer)
    zout.close()
    zin.close()

    try:
        print( subprocess.check_output(["cmd", "/c", "copy", zip_path, backup_dir]) )
        zf = zipfile.ZipFile(os.path.join(backup_dir, os.path.join(backup_dir,filename)))
        print ("Backup to external drive successful!")
    except:
        print ("Backup to external drive did not work!")
        raise

