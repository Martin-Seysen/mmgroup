import sys
import re

sys.path.append(r".")
from find_generators import ORDER7_FILENAME










cases_pattern = r"\s*\#\s*([0-9]+)\s+cases"
cases_match = re.compile(cases_pattern)

time_pattern = r"\s*\#\s*Time:\s*[-A-Za-z0-9.:]+\s+,\s+duration\s+([0-9.]+)\s+seconds"
time_match = re.compile(time_pattern)




def parse_datafile(filename, verbose = 0):
    print("Reading data from file '%s'..." % filename)
    n_cases = 0
    time_ = 0.0
    f = open(filename)
    for line in f.readlines():
         m = cases_match.match(line)
         if m:
             if verbose:
                 print(m.groups())
             if time1 is not None:
                 n_cases += int(m.groups()[0])
                 time_ += time1
         m = time_match.match(line)
         if m:
             if verbose:
                 print(m.groups())
             time1 = float(m.groups()[0])
         else:
             time1 = None
    if verbose:
        print("time: %.1f s, %d cases" % (time_, n_cases ))
        #print(data)
    return time_, n_cases



if __name__ == "__main__":
    FILENAME = ORDER7_FILENAME
    time_, n_cases = parse_datafile(FILENAME)
    MegaCasesPerDay = (n_cases / time_ ) * (86400 / 1.0e6) * 0.95
    s = """The computer that has created file '%s'
can check about %.1f million cases per day."""
    print(s % (ORDER7_FILENAME, MegaCasesPerDay))


