"""Create a temporary directory for tests

For some tests we want to create a temporary directory in a
known position in the source file tree, so that C programs or 
python extensions can be created end executed in that temporary 
directory.

Function make_tmp_dir(parent) creates such a temporary directory 
and returns its path name. Here 'parent' can be a path name or a 
python module object. Function kill_tmp_dir(parent) deletes that 
directory. Function temp_dir_name(parent) returns the path name
of the temporary directory without creating it. 

If '<parent>' is a python module object, then '<parent>.temp' 
is a module object referring to the temporary directory created. 
"""

import os
import shutil
import time
import types

TMPNAME = "temp"

temp_dir = None

def temp_dir_name(parent= None):
    """Return absolute path name of temporary directory

    Here argument 'parent' is as in function make_temp_dir().
    This function does not create any directory.
    """
    if isinstance(parent, types.ModuleType):
        dir_ = os.path.dirname(os.path.realpath(parent.__file__))
        return  os.path.abspath(os.path.join(dir_, TMPNAME))
    elif parent:
        return os.path.abspath(os.path.join(parent, TMPNAME))
    else:
        return os.path.abspath(TMPNAME)


def make_temp_dir(parent = None): 
    """Create a temporary directory, return its absolute path name

    Optionally, a 'parent' directory of the temporary directory
    may be given. By default, this is the current directory.

    Alternatively, 'parent' be a module object which has been 
    imported with an "import" statement. Then a directory 'temp' 
    is created as a subdirectory of the directory of that module, 
    which is os.path.dirname(os.path.realpath(module.__file__)).
    Futhermore, and empty file with name "__init__.py" is created
    in that subdirectory 'temp'. Then a python module foo.py
    written into that subdirectory can be imported as
    <parent>.temp.foo, with <parent> the module object 'parent'.

    Caution: this destroys an existing file or directory with name
    "temp" in the parent directory.
    """ 
    global temp_dir  
    temp_dir =  temp_dir_name(parent)
    try:
        shutil.rmtree(temp_dir)
        time.sleep(0.1) # It take a while to delete TMPNAME
    except FileNotFoundError:
        pass
    os.mkdir(temp_dir)
    if isinstance(parent, types.ModuleType):
        f = open(os.path.join(temp_dir, "__init__.py"), "wt")
        f.write("\n")
        f.close()
    return temp_dir


def kill_temp_dir(name):
    """Remove a temporary directory

    Here 'name' must be the path name returned from the last
    recent call to function make_tmp_dir().
    """
    global temp_dir  
    if name and name == temp_dir:
        try:
            shutil.rmtree(temp_dir)
        except FileNotFoundError:
            pass
        temp_dir = None
    


    