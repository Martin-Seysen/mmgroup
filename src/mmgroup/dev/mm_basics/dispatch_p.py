import sys
import os
import re
from collections import defaultdict

from mmgroup.generate_c import UserDirective, EmptyUserDirective
from mmgroup.generate_c import UserFormat, ZeroUserFormat

sys.path.append(os.path.join('..','..','..','..'))
from config import C_DIR
sys.path.pop()


m_files = re.compile("mm_op(\d+)\.h")

def iter_find_files(dir):
    for fname in os.listdir(dir):
        m = m_files.match(fname)
        if m:
            yield int(m.group(1)), os.path.join(dir, fname)


m_declaration = re.compile(
  "\s*(int\d+_t|void)\s+(mm_op([0-9]+)[A-Za-z0-9_]*)\s*\("
)


def function_names(c_dir):
    d = {}
    for p, file_name in  iter_find_files(c_dir):
        d[p] = {}
        for line in open(file_name, "rt"):
            m = m_declaration.match(line)
            if (m):
                return_type =  m.group(1)
                function_name = m.group(2)
                p_in_function_name = int(m.group(3))
                if p == p_in_function_name:
                    d[p][function_name] = return_type
    return d



def make_p_case(p, function, return_type, *args):
    s = "        case %d:\n" % p
    invocation = "%s(%s);\n" % (function, ", ".join(args))
    if return_type == "void":
        s += "            " + invocation;
        s += "            return 0;\n"
    else:
        s += "            return " + invocation
    return s


def make_C_switch(primes, return_type, function_name, p_var, *args):
    if function_name[:5] != "mm_op":
         raise ValueError("Function name must satrt with mm_op")

    s = "    switch (%s) {\n" % p_var
    for p in primes:
        called_function =  function_name[:5] + str(p) + function_name[5:]
        #print(called_function)
        s +=  make_p_case(p, return_type, called_function, *args)
    s += """        default:
            return -1;
    }
"""
    return s



m_files = re.compile("(mm_op)(\d+)(.+)")


class DispatchP:
    def  __init__(self, **kwds):
        c_files = getattr(kwds, 'c_files', "")
        if isinstance(c_files, list):
             c_files = "\n".join(c_files)
        d = defaultdict(list)
        for name in c_files.split():
            m = mfiles.match(name)
            if m:
                 d[m.groups(0) + m.groups(2)].append(m.groups(1))
        self._d = dict(d)
   

    def dispatch_p(self, arg_string):
        args = [s.strip() for s in arg_string.split(',')] 
        return make_C_switch(self.d, *args)

    def legal_p(self, function_name):
        if function_name[:5] != "mm_op":
            raise ValueError("Function name must start with 'mm_op'")

        p_list = sorted(getattr(self._d, function_name, [])[:])
        if len(p_list > 1):
             p_list[-1] = "or " + p_list[-1] 
        return ", ".join(str(p) for p in sorted(p_list))
       
    @property
    def tables(self):
        return {
            'LEGAL_P': UserFormat(self.legal_p, 's')
        }
       
    @property
    def directives(self):
        return {
            'DISPATCH_P': UserDirective(self.dispatch_p, '.')
        }
      
    mockup_directives = {'DISPATCH_P': EmptyUserDirective}


Tables = DispatchP


#print( make_C_switch(D, 'mm_op_pi', 'p', 'a', 'b', 'c') )
#print( make_C_switch(D, 'mm_op_load_leech3matrix', 'p', 'a', 'b', 'c') )


