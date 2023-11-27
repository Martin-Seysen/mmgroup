import sys
import os
import re
from collections import defaultdict

sys.path.append(os.path.join('..','..','..'))
sys.path.append(os.path.join('..','..','..','..'))

from mmgroup.generate_c import UserDirective, EmptyUserDirective
from mmgroup.generate_c import UserFormat, ZeroUserFormat
from mmgroup.generate_c import UserFormat, ZeroUserFormat
from mmgroup.generate_c import iter_exports_from_header





#####################################################################
# Extract declarations from .h files
#####################################################################

#H_FILES = ["mm_op%d.h" % p for p in (3, 7, 15, 31, 127, 255)]

H_FILES = ["mm_op_sub.h"]


m_declaration = re.compile(
  r"\s*(int\d+_t|void)\s+(mm_op([0-9]+)[A-Za-z0-9_]*)\s*\("
)

def iter_find_export(c_dir):
    """Scan prototypes from headers ``H_FILE`` in directory ``c_dir``

    The function yields triples of strings
    ``(return_type, function_name, modulus)``
    of relevant functions scanned from the header files
    """
    for h_name in H_FILES:
        exp = iter_exports_from_header(os.path.join(c_dir, h_name))
        for kwd, args, prototype in exp:
            prototype = prototype.strip()
            if kwd == "EXPORT" :
                m = m_declaration.match(prototype)
                if m:
                    ret_type, name, p = m.groups()
                    raw_name = re.sub(r"mm_op([0-9]+)", "mm_op", name)
                    yield ret_type, raw_name, p
       




def function_names(c_dir):
    d = {}
    for return_type, function_name, p in  iter_find_export(c_dir):
        if not function_name in d:
            d[function_name] = (return_type, []) 
        assert d[function_name][0] == return_type
        d[function_name][1].append(int(p))
    for key in d:
        d[key][1].sort()
    return d



#####################################################################
#  generate C switch / case statements
#####################################################################


def make_p_case(p, return_type, function, *args):
    s = "        case %d:\n" % p
    invocation = "%s(%s);\n" % (function, ", ".join(args))
    if return_type == "void":
        s += "            " + invocation;
        s += "            return 0;\n"
    else:
        s += "            return " + invocation
    return s


def make_C_switch(c_functions, return_type, function_name, p_var, *args):
    if function_name[:5] != "mm_op":
         raise ValueError("Function name must satrt with mm_op")
    expected_return_type, primes = c_functions[function_name]
    assert return_type == expected_return_type
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

#####################################################################
#  Creating tables an directives
#####################################################################




class DispatchP:
    c_functions = None
    def  __init__(self, **kwds):
        c_dir = kwds.get('C_DIR', None)
        if c_dir:
            self.c_functions = function_names(c_dir)
   

    def dispatch_p(self, arg_string):
        args = [s.strip() for s in arg_string.split(',')] 
        return make_C_switch(self.c_functions, *args)

    def legal_p_text(self, function_name):
        if function_name[:6] != "mm_op_":
            raise ValueError("Function name must start with 'mm_op_'")
        p_list = [str(p) for p in self.c_functions[function_name][1]]
        if len(p_list) > 1:
            p_list[-1] = 'or ' + p_list[-1] 
        return ", ".join(p_list)

    def p_table(self):
        pl = self.c_functions['mm_op_word'][1] if self.c_functions else [] 
        return pl + [0]
       
    @property
    def tables(self):
        return {
            'LEGAL_P': UserFormat(self.legal_p_text, 's'),
            'P_LIST': self.p_table(),
        }
       
    @property
    def directives(self):
        return {
            'DISPATCH_P': UserDirective(self.dispatch_p, '.')
        }
      


Tables = DispatchP


#################################################################
# Creating legacy python modules mm_<p>.py
#################################################################


# List of existing C functions obtained by excuting
# function iter_find_export('...\src\mmgroup\dev\c_files')
# at 2023-07-10. This will be used for creating legacy code.
ALL_LEGACY_FUNCTIONS = {
'mm_op_pi': ('int32_t', [3, 7, 15, 31, 127, 255]), 
'mm_op_copy': ('int32_t', [3, 7, 15, 31, 127, 255]), 
'mm_op_compare_len': ('int32_t', [3, 7, 15, 31, 127, 255]), 
'mm_op_compare': ('int32_t', [3, 7, 15, 31, 127, 255]), 
'mm_op_checkzero': ('int32_t', [3, 7, 15, 31, 127, 255]), 
'mm_op_vector_add': ('int32_t', [3, 7, 15, 31, 127, 255]), 
'mm_op_scalar_mul': ('int32_t', [3, 7, 15, 31, 127, 255]), 
'mm_op_compare_mod_q': ('int32_t', [3, 7, 15, 31, 127, 255]), 
'mm_op_store_axis': ('int32_t', [3, 7, 15, 31, 127, 255]), 
'mm_op_xy': ('int32_t', [3, 7, 15, 31, 127, 255]), 
'mm_op_omega': ('int32_t', [3, 7, 15, 31, 127, 255]), 
'mm_op_t_A': ('int32_t', [3, 7, 15, 31, 127, 255]), 
'mm_op_word': ('int32_t', [3, 7, 15, 31, 127, 255]), 
'mm_op_word_tag_A': ('int32_t', [3, 7, 15, 31, 127, 255]), 
'mm_op_word_ABC': ('int32_t', [3, 7, 15, 31, 127, 255]),
'mm_op_eval_A_rank_mod3': ('int64_t', [3, 15]), 
'mm_op_load_leech3matrix': ('int32_t', [3, 15]), 
'mm_op_eval_A_aux': ('int32_t', [3, 15]), 
'mm_op_eval_A': ('int32_t', [3, 15]), 
'mm_op_norm_A': ('int32_t', [3, 15]), 
'mm_op_watermark_A': ('int32_t', [3, 15]),
'mm_op_watermark_A_perm_num': ('int32_t', [3, 15]), 
'mm_op_eval_X_find_abs': ('int32_t', [15]), 
'mm_op_eval_X_count_abs': ('int32_t', [15])
}

ALL_PRIMES = [3, 7, 15, 31, 127, 255]


def write_f_mapping(f, name, p):
    s = """def {name}(*args, **kwds):
    return mm_op.mm_{name}({p}, *args, **kwds)

"""
    f.write(s.format(name = name, p = p))
  
     


def make_one_legacy_script(out_dir, p):
    f_name = os.path.join(out_dir, "mm%d.py" % p)
    f = open(f_name, "wt")
    s = '''"""This module is deprecated; do not use in new projects!"""

import warnings
from mmgroup import mm_op

warnings.warn("Module mmgroup.mm{0} is deprecated! " 
    "Replace function 'op_<xxx>(*args)' in this module by function "
    "'mm_op_<xxx>({0}, *args)' in module mmgroup.mm_op!",
    UserWarning)

'''

    f.write(s.format(p))

    for name, (return_type, primes) in ALL_LEGACY_FUNCTIONS.items():
        if p in primes:
            assert name[:3] == 'mm_'
            write_f_mapping(f, name[3:], p)
    f.close() 

     
def make_legacy_scripts_mod_p(out_dir):
    for p in ALL_PRIMES:
        make_one_legacy_script(out_dir, p)


#################################################################
# Creating legacy python module mm.py
#################################################################


MM_VARS = [
 'PROTECT_OVERFLOW', 'INT_BITS',
 'mm_aux_mmv_size',
 'mm_aux_zero_mmv',
 'mm_aux_random_mmv', 'mm_aux_reduce_mmv', 'mm_aux_reduce_mmv_fields',
 'mm_aux_check_mmv', 'mm_aux_small24_expand', 'mm_aux_small24_compress',
 'mm_aux_mmv_to_bytes', 'mm_aux_bytes_to_mmv', 'mm_aux_mmv_to_sparse', 
 'mm_aux_mmv_extract_sparse', 'mm_aux_mmv_get_sparse',
 'mm_aux_mmv_add_sparse', 'mm_aux_mmv_set_sparse',
 'mm_aux_mmv_extract_sparse_signs', 'mm_aux_mmv_extract_x_signs',
 'mm_aux_mul_sparse', 'mm_aux_index_extern_to_sparse',
 'mm_aux_array_extern_to_sparse', 'mm_aux_index_sparse_to_extern',
 'mm_aux_index_sparse_to_leech', 'mm_aux_index_sparse_to_leech2',
 'mm_aux_index_leech2_to_sparse', 'mm_aux_index_intern_to_sparse',
 'mm_aux_hash', 'mm_sub_test_prep_pi_64',
 'mm_sub_test_prep_xy', 'mm_group_prepare_op_ABC',
 'mm_sub_get_table_xi', 'mm_crt_combine', 'mm_crt_combine_bytes',
 'mm_crt_check_v2', 'mm_crt_check_g', 'mm_crt_norm_int32_32',
 'mm_crt_norm_int32', 'mm_crt_v2_int32_32', 'mm_crt_v2_int32',
 'mmv_array', 'mm_vector'
]


def make_mm_legacy_script(out_dir):
    f_name = os.path.join(out_dir, "mm.py")
    print(f_name)
    f = open(f_name, "wt")
    s = '''"""This module is deprecated; do not use in new projects!"""

import warnings
from mmgroup import mm_op

warnings.warn("Module mmgroup.mm is deprecated! " 
    "Replace a function in this module by the "
    "corresponding function in module mmgroup.mm_op!",
    UserWarning)


def mm_aux_read_direct_mmv1(*args, **kwds):
    return mm_op.mm_aux_read_mmv1(*args, **kwds) 

'''
    f.write(s)
    for name in MM_VARS:
        f.write("{0} = mm_op.{0}\n".format(name))
    f.close() 


#################################################################
# Main program
#################################################################



def make_all_legacy_scripts(out_dir):
    make_mm_legacy_script(out_dir)
    make_legacy_scripts_mod_p(out_dir)




if __name__ == "__main__":
    make_all_legacy_scripts(sys.argv[1])

