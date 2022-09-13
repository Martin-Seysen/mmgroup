import sys
import re
import argparse

sys.path.append(r".")
from find_generators import ORDER7_FILENAME
from check import HurwitzVerifyer

HURWITZ_MONSTER_FILENAME = "hurwitz_monster_samples.py"








mm_pattern = r"\s*(\w+)\s*=\s*\'?(M0\<[A-Za-z0-9_* ]+\>)\'?"
mm_match = re.compile(mm_pattern)

order_pattern = r"\s*order\s*=\s*([0-9]+)"
order_match = re.compile(order_pattern)

cases_pattern = r"\s*\#\s*([0-9]+)\s+cases"
cases_match = re.compile(cases_pattern)


def parse_order7_datafile(filename, verbose = 0):
    print("Reading data from file %s..." % filename)
    g3 = ge = None
    data = []
    n_cases = 0
    f = open(filename)
    for line in f.readlines():
         m = mm_match.match(line)
         if m:
             if verbose:
                 print(m.groups())
             if m.groups()[0] == 'g3': g3 = m.groups()[1]
             if m.groups()[0] == 'ge': ge = m.groups()[1]
         m = order_match.match(line)
         if m:
             #print(m.groups())
             order = int(m.groups()[0])
             if order == 7 and g3 is not None and ge is not None:
                 data.append((g3, ge))
             g3 = ge = None
         m = cases_match.match(line)
         if m:
             n_cases += int(m.groups()[0])
             g3 = ge = None
    if verbose:
        print(n_cases, "cases")
        #print(data)
    return data, n_cases



HURWITZ_MONSTER_FILENAME_HEADER = (
"""# This file has been generated automatically, do not change!
# It contains sets of generators of the Monster as a Huwitz group.
# See project documentation for the interpretation of these data.

hurwitz_monster_data = []
"""
)

N_CASES = r"""
# Generators of the Monster as a Hurwitz group 
# found after %d trials
"""

CASE_STR = r"""
g3 ='%s'
ge ='%s'
verifiers = [0x%016x, 0x%016x]
hurwitz_monster_data.append((g3, ge, verifiers))
"""


def write_Hurwitz_Monster_datafile(data, n_cases, verbose = 0):
    print("Looking for a Hurwitz group presentation of the Monster")
    data_string_list = []
    for g3, ge in data:
        if verbose:
            print("Checking input sample...")
        ver = HurwitzVerifyer(g3, ge)
        if ver.verify(n_trials = 400, verbose = verbose):
            verifiers = list(ver.verifiers.values())
            if len(verifiers) == 2:
                #print(g3, ge, verifiers[0], verifiers[1])
                s = CASE_STR % (g3, ge, verifiers[0], verifiers[1])
                data_string_list.append(s)
                if verbose:
                    print("Hurwitz sample generates the Monster")
            else:
                print("Bad Hurwitz group sample")
        else:
            if verbose:
                print("Hurwitz sample does not generate the Monster")

 
    print("Writing output to file " + HURWITZ_MONSTER_FILENAME)
    try:
        f = open(HURWITZ_MONSTER_FILENAME, "rt")
        f.close()
        f = open(HURWITZ_MONSTER_FILENAME, "at")
    except:
        f = open(HURWITZ_MONSTER_FILENAME, "at")
        f.write(HURWITZ_MONSTER_FILENAME_HEADER)
    f.write(N_CASES % n_cases)
    if (len(data_string_list)):
        f.write("\n".join(data_string_list) + "\n")
    else:
        f.write("# (No data found)\n\n")
    f.close()
       


def make_parser():
    description = """Read output of the script find_generators.py and write 
sets of generators of the Monster as a Hurwitz group to file %s.
"""

    parser = argparse.ArgumentParser(description
        = description % HURWITZ_MONSTER_FILENAME )
    help_input = "Name of input file, default is '%s'"
    parser.add_argument('InputFile',
                       metavar='InputFile',
                       type=str,
                       default=ORDER7_FILENAME,
                       nargs='?',
                       help=help_input % ORDER7_FILENAME )
    parser.add_argument("-v",  dest="verbose", action="store_true",
        help="Verbose operation" )
    return parser


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args()
    verbose = args.verbose
    in_file =  args.InputFile
    #print(verbose, in_file)
    data, n_cases = parse_order7_datafile(in_file, verbose = 0)
    write_Hurwitz_Monster_datafile(data, n_cases, verbose = verbose)


