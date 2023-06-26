import sys
import os
import re

from mmgroup.generate_c import UserDirective, UserFormat

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


D = function_names(C_DIR)

def make_C_switch(d, function_name, p_var, *args):
    if function_name[:5] != "mm_op":
         raise ValueError("Function name must satrt with mm_op")

    s = "    switch (%s) {\n" % p_var
    for p, f_dict in sorted(d.items()):
        called_function =  function_name[:5] + str(p) + function_name[5:]
        #print(called_function)
        if called_function  in f_dict:
            return_type = f_dict[called_function]
            s +=  make_p_case(p, called_function, return_type, *args)
    s += """        default:
            return -1;
    }
"""
    return s





class DispatchP():
    def  __init__(self, c_dir):
        self.c_dir = c_dir
        self._d = None
   
    @property
    def d(self):
        if self._d is None:
            self._d = function_names(self.c_dir)
        return self._d


    def dispatch_p(self, arg_string):
        args = [s.strip() for s in arg_string.split(',')] 
        return make_C_switch(self.d, *args)

    def legal_p(self, function_name):
        if function_name[:5] != "mm_op":
            raise ValueError("Function name must satrt with mm_op")

        p_list = []
        for p, f_dict in self.d.items():
            called_function =  function_name[:5] + str(p) + function_name[5:]
            if called_function in f_dict:
                p_list.append(p)
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
      


class  Mockup_DispatchP(DispatchP):
    def __init__(self, c_dir):
        super(Mockup_DispatchP, self).__init__(c_dir)

    def dispatch_p(self, *args):
        return ""



#print( make_C_switch(D, 'mm_op_pi', 'p', 'a', 'b', 'c') )
#print( make_C_switch(D, 'mm_op_load_leech3matrix', 'p', 'a', 'b', 'c') )


