"""Generate code for computing scalar product in the rep of the monster

"""

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

from mmgroup.generate_c import UserDirective, EmptyUserDirective

###########################################################################
# Computing scalar product (mod p) of vector in the rep of the monster
###########################################################################


def gen_scalprod_mod3_snippet(is24):
    """Scalar multiplication mod 3

    Return a code snipped the computes "_ac1", "_ac2" with
    2 * "_ac2" + "_ac1" = scalar_product("_v1", "_v2") where "_v1",
    "_v2" are 64-bit integers containing a vector of 32 integers mod 3.
    Up to 0x7fff such values can be summed up without intermediate
    reduction. If is24 is True then "_v1", "_v2" are 64-bit integers
    containing a vector of 24 integers mod 3. The function uses a
    temporary 64-bit integer "_t". The names of the variables cannot
    be changed.
    """
    if  is24:
        lowmask = 0x555555555555
        tvalue =  '_t'
    else:  
        lowmask =  0x5555555555555555
        tvalue = '(_t >> 32) + (_t & 0xffffffffULL)'
    highmask = 2*lowmask
    return f"""    _t = _v1 & 0x{lowmask:x}ULL;
    _t = (_t + _t + _t) & _v2;
    _ac1 += {tvalue};
    _t = _v1 & 0x{highmask:x}ULL;
    _t = (_t + (_t >> 1)) & _v2;
    _ac2 += {tvalue};"""


def gen_scalprod_mod3(mv1, mv2, length, is24):
    scalprod = gen_scalprod_mod3_snippet(is24)
    return f"""uint_mmv_t _v1, _v2, _t, _ac1 = 0, _ac2 = 0;
do {{
    _v1 = *{mv1}++;
    _v2 = *{mv2}++;
{scalprod}
}} while (--{length});
return (uint32_t)((_ac2 + _ac2 + _ac1) % 3);
"""


def gen_scalprod_mod3_index(mv1, mv2, ind, imax):
    scalprod = gen_scalprod_mod3_snippet(False)
    return f"""uint_mmv_t _v1, _v2, _t, _ac1 = 0, _ac2 = 0;
uint32_t _i;
while(1) {{
    _i = *((*{ind})++);
    if (_i >= {imax}) break;
    _v1 = {mv1}[_i];
    _v2 = {mv2}[_i];
{scalprod}
}}
--(*{ind});
return (uint32_t)((_ac2 + _ac2 + _ac1) % 3);
""" 



def gen_extract(v1, v2, res, width, bit):
    mask = 0x1111111111111111 << bit
    shl = width - bit
    shr = f"({res} >> {bit})" if bit else res
    C1 = f"extract entry of {v2} if bit {bit} of entry in {v1} is set"
    return f"""    // {C1} 
    {res} = {v1} & 0x{mask:x}ULL;
    {res} = ({res} << {shl}) - {shr};
    {res} &= {v2}; // extracted entries of {v2} in {res}
"""


def gen_adddown(v1, width, cy=False):
    """Put v1 =  v1_odd + v1_even

    Here the 64-bit integer v1 is considered as an array of integers
    of with ``width``, indexed from 0 to 64 / width -1.
    The function adds the entries with index 2*i and 2*i + 1 and
    stores the result at entry 2*i, with the sum possibly one bit
    longer. If ``cy`` is set to 0 then some masking is saved; details
    are hairy! If in doubt, set ``cy = 1``.
    """
    if width == 32:
       return f"    {v1} = ({v1} >> 32) + ({v1} & 0xffffffffULL);\n"
    m = 0xffffffffffffffff // ((1 << (2*width)) - 1)
    m *= (1 << width) - 1
    if cy:
       s = f"""({v1} & 0x{m:x}ULL) 
      + (({v1} >> {width}) & 0x{m:x}ULL)""" 
    else:
       s =  f"({v1} + ({v1} >> {width})) & 0x{m:x}ULL"
    return f"    {v1} = {s};\n"


def gen_add_scal_mod7_snippet(v1, v2, res):
    return "".join([
       gen_extract(v1, v2, "_t1", 3, 2),
       gen_adddown("_t1", 4),
       gen_extract(v1, v2, "_t2", 3, 1),
       gen_adddown("_t2", 4),
       "    _t1 = _t1 + _t1 + _t2;\n",
       gen_extract(v1, v2, "_t2", 3, 0),
       gen_adddown("_t2", 4),
       "    _t1 = _t1 + _t1 + _t2;\n",
       gen_adddown('_t1', 8),
       f"    {res} += _t1;\n",
    ])


def gen_scalprod_mod7(mv1, mv2, length):
    return "".join([ f"""uint_mmv_t _v1, _v2, _res0, _res1 = 0, _t1, _t2;
uint32_t _i, _j;
{length} <<= 1;
_i = ((uint32_t)length + 0xffULL) & ~0xffULL;
_j = 0x100ULL + (uint32_t)length  - _i;
_i >>= 8;  // Here {length} = 0x100 * (_i - 1) + _j; _j > 0 
do {{
    _res0 = 0;
  do {{
    _v1 = *({mv1}++);
    _v2 = *({mv2}++);
""",
    gen_add_scal_mod7_snippet('_v1', '_v2', '_res0'),
    "  } while (--_j);\n",
    gen_adddown('_res0', 16, cy=True),
    gen_adddown('_res0', 32),
    "    _res1 += _res0;\n",
    "    _j = 0x100;\n",
    "} while (--_i);\n",
    "return (uint32_t)(_res1 % 7);\n"
]) 
      

def gen_add_scal_mod15_snippet(v1, v2, res):
    return "".join([
       gen_extract(v1, v2, "_t1", 4, 3),
       gen_adddown("_t1", 32),
       gen_extract(v1, v2, "_t2", 4, 2),
       gen_adddown("_t2", 32),
       "    _t1 = _t1 + _t1 + _t2;\n",
       gen_extract(v1, v2, "_t2", 4, 1),
       gen_adddown("_t2", 32),
       "    _t1 = _t1 + _t1 + _t2;\n",
       gen_extract(v1, v2, "_t2", 4, 0),
       gen_adddown("_t2", 32),
       "    _t1 = _t1 + _t1 + _t2;\n",
       f"    {res} += _t1;\n",
    ])


def gen_scalprod_mod15(mv1, mv2, length):
    return "".join([ f"""uint_mmv_t _v1, _v2, _res = 0, _t1, _t2;
{length} <<= 1;
do {{
    _v1 = *({mv1}++);
    _v2 = *({mv2}++);
""",
    gen_add_scal_mod15_snippet('_v1', '_v2', '_res'),
    f"}} while (--{length});\n",
    "return (uint32_t)(_res % 15);\n"
]) 


def gen_scalprod_mod15_index(mv1, mv2, ind, imax):
    scalprod = gen_add_scal_mod15_snippet('_v1', '_v2', '_res')
    return f"""uint_mmv_t _v1, _v2, _res = 0, _t1, _t2;
uint32_t _i;
while(1) {{
    _i = *((*{ind})++);
    if (_i >= {imax}) break;
    _v1 = {mv1}[_i];
    _v2 = {mv2}[_i];
{scalprod} }}
--(*{ind});
return (uint32_t)(_res % 15);
"""





def gen_add_scal_mod_large_snippet(pbits, v1, v2, res):
    assert 1 <= pbits <= 8
    if pbits == 8:
        s = [f"    {res} += ({v1} & 0xff) * ({v2} & 0xff)\n"]
        for i in range(8, 64, 8):
            t = ";\n" if i == 56 else "\n"
            s.append(
            f"    + (({v1} >> {i}) & 0xff) * (({v2} >> {i}) & 0xff){t}"
            )
        return "".join(s)
    p = (1 << pbits) - 1
    p2 = (1 << (2*pbits)) - 1
    pd = hex(p + (p << 16)) + 'ULL'
    p2d = hex(p2 + (p2 << 32)) + 'ULL'
    return f"""    {res} += (({v1} & {pd}) * ({v2} & {pd})
      & {p2d}) +
    ((({v1} >> 8) & {pd}) * (({v2} >> 8) & {pd})
      & {p2d}) +
    ((({v1} >> 32) & {pd}) * (({v2} >> 32) & {pd})
      & {p2d}) +
    ((({v1} >> 40) & {pd}) * (({v2} >> 40) & {pd})
      & {p2d});
"""


def gen_add_scal_mod_large_32(pbits, mv1, mv2, length):
    p = (1 << pbits) - 1
    s = [ f"""uint_mmv_t _v1, _v2, _res = 0;
{length} <<= 2;
do {{
    _v1 = *({mv1}++);
    _v2 = *({mv2}++);
""" ,
    gen_add_scal_mod_large_snippet(pbits, '_v1', '_v2', '_res'),
    f"""}} while (--{length});\n"""
    ]
    if pbits < 8:
        s.append(gen_adddown("_res", 32))
    s.append(f"return (uint32_t)(_res % {p});\n")
    return  "".join(s)



def gen_add_scal_mod_large_24(pbits, mv1, mv2, length):
    p = (1 << pbits) - 1
    s = [ f"""uint_mmv_t _v1, _v2, _res = 0, _i;
do {{
    _i = 3;
  do {{
    _v1 = *({mv1}++);
    _v2 = *({mv2}++);
""" ,
    gen_add_scal_mod_large_snippet(pbits, '_v1', '_v2', '_res'),
    f"""  }} while (--_i);
    ++{mv1}; ++{mv2};
}} while (--{length});
"""
    ]
    if pbits < 8:
        s.append(gen_adddown("_res", 32))
    s.append(f"return (uint32_t)(_res % {p});\n")
    return  "".join(s)


######################################################################
# Table classes
######################################################################



class Tables: #(MM_Op):
    BIT_DICT = {3:2, 7:3, 15:4, 31:5, 63:6, 127:7, 255:8}
    def __init__(self, **kwds):
        #super(Tables, self).__init__(**kwds)
        self.p = p = int(kwds.get('p', 3))
        assert self.p in [3, 7, 15, 31, 63, 127, 255]
        self.pbits = self.BIT_DICT[self.p]
        self.tables = {
            "HAS_SCALPROD_MOD_P_24" : not self.p in [7, 15]
        }
        self.directives = {
            "SCALPROD_MOD_P_32" : UserDirective(self.scalprod_32, "sss"),
            "SCALPROD_MOD_P_24" : UserDirective(self.scalprod_24, "sss"),
            "SCALPROD_MOD_P_INDEX" :
                              UserDirective(self.scalprod_index, "ssss"),
        }

    def scalprod_32(self, mv1, mv2, length):
        if self.p == 3:
            return gen_scalprod_mod3(mv1, mv2, length, False)
        elif self.p == 7:
            return gen_scalprod_mod7(mv1, mv2, length)
        elif self.p == 15:
            return gen_scalprod_mod15(mv1, mv2, length)
        else:
            return gen_add_scal_mod_large_32(self.pbits, mv1, mv2, length)

    def scalprod_24(self, mv1, mv2, length):
        if self.p == 3:
            return gen_scalprod_mod3(mv1, mv2, length, True)
        elif self.p == 7:
            return gen_scalprod_mod7(mv1, mv2, length)
        elif self.p == 15:
            return gen_scalprod_mod15(mv1, mv2, length)
        else:
            return gen_add_scal_mod_large_24(self.pbits, mv1, mv2, length)

    def scalprod_index(self, mv1, mv2, ind, imax):
        if self.p == 3:
            return gen_scalprod_mod3_index(mv1, mv2, ind, imax)
        elif self.p == 15:
            return gen_scalprod_mod15_index(mv1, mv2, ind, imax)
        ERR = "Modulus %d not supported for directive SCALPROD_MOD_P_INDEX"
        raise ValueError(ERR % self.p)



class MockupTables:
    tables = {
        "HAS_SCALPROD_MOD_P_24" : False
    }
    directives = {
        "SCALPROD_MOD_P_32" : EmptyUserDirective,
        "SCALPROD_MOD_P_24" : EmptyUserDirective,
        "SCALPROD_MOD_P_INDEX" : EmptyUserDirective,
    }
    def __init__(self, **kwds):
        pass

######################################################################
# Test functions
######################################################################

if __name__ == "__main__":
    1/0





