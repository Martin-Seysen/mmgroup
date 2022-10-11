"""Compute Chisquare statistics

The main function exported from this module is function chisquare, which
is a substitute for function scipy.special.gammainc in the scipy package.

Rewriting a fast implementation in slow python code appears to be stupid,
but downloading about ten versions of the scipy package for building ten
different wheels is even slower.

For a regression test of the functions in this module against the
corresponding functions in the scipy package, one should execute this 
python module in a shell.

Caution:

Although there is hardly any need for a change, the regression test 
must be run manually after any change in this module!


For background see:

[1] William H. Press and Saul A. Teukolsky and William T. Vetterling 
    and Brian P. Flannery, Numerical Recipes in C, Cambridge University 
    Press, 1992.

[2] W. Gautschi, A Computational Procedure for the inconplete Gamma Functions
    https://www.cs.purdue.edu/homes/wxg/selected_works/section_02/068.pdf



"""

import math
from math import log, exp, fabs, lgamma

ERR_INC_GAMMA_RANGE = "Invalid arguments in incomplete gamma funtion"
ERR_INC_GAMMA_A_LARGE = "a is too large in incomplete gamma funtion"

ITMAX = 100
EPS = 3.0e-7
FPMIN = 1.0e-30

# For sample code of low-level functions see: 
# https://personal.sron.nl/~kuiper/hea_inst11/exercise-2/Varia/Numerical-Recipes/f6-2.pdf

def gamma_ser(a, x):
    """Incomplete Gamma function computes with power series

    Not for public use!
    This is equivalent to scipy.special.gammainc(a, x).
    It is efficient for small values x, see [1], section 6.2.
    """   
    a, x = a+0.0, x+0.0
    if x <= 0.0 or a <= 0.0:
        if x == 0.0 and a > 0.0: return 0.0
        raise ValueError(ERR_INC_GAMMA_RANGE)
    gln = lgamma(a)
    ap = a
    summ = 1.0 / a
    dele = summ
    for n in range(ITMAX):
        ap += 1.0
        dele *= x/ap
        summ += dele
        if fabs(dele) < fabs(summ)*EPS:
            return summ * exp(-x + a *log(x) - gln)
    raise ValueError(ERR_INC_GAMMA_A_LARGE)



def gamma_cfrac(a, x):
    """Incomplete Gamma function computes with continued fractions

    Not for public use!
    This is equivalent to scipy.special.gammainc(a, x).
    It is efficient for large values x, see [1], section 6.2.
    """
    a, x = a+0.0, x+0.0
    if x <= 0.0 or a <= 0.0:
        raise ValueError(ERR_INC_GAMMA_RANGE)
    gln = lgamma(a)
    b = x + 1.0 - a
    c = 1.0/FPMIN
    d = 1.0/b
    h = d
    for i in range(1, ITMAX+1):
        an = -i * (i-a)
        b += 2.0
        d = an*d + b
        if (fabs(d) < FPMIN): d = FPMIN
        c  = b + an / c
        if (fabs(c) < FPMIN): c = FPMIN
        d = 1.0/d
        dele = d*c
        h *= dele
        if (fabs(dele-1.0) < EPS):
            return exp(-x + a*log(x) - gln) * h
    raise ValueError(ERR_INC_GAMMA_A_LARGE)


def gamma_p(a, x):
    """Incomplete Gamma function 

    This is equivalent to scipy.special.gammainc(a, x).
    """
    if x < (a + 1.0):
         return gamma_ser(a, x)
    else:
         return 1.0 - gamma_cfrac(a, x)


def gamma_q(a, x):
    """Complement of incomplete Gamma function 

    This is equivalent to  1.0 - scipy.special.gammainc(a, x)
    """
    if x < (a + 1.0):
         return 1.0 - gamma_ser(a, x)
    else:
         return gamma_cfrac(a, x)
     

def chisquare_p(chisq, nu):
    """The (complemented) chisquare function see [1], section 6.2"""
    return gamma_q(0.5 * nu, 0.5 * chisq)


def chisquare(f_obt, f_exp = None):
    """Chisquare test

    The lists ``f_obt`` and ``f_exp`` describe the obtained
    and the expected distribution. ``f_exp`` defaults to the
    uniform distibution.

    The function returns a pair ``(chisq, p)``. This function
    is equivalent to function ``scipy.stats.chisquare()``.
    """
    n = len(f_obt)
    if f_exp is None:
        f_exp = [1.0] * n
    assert len(f_exp) == n > 1
    assert min(f_exp) > 0
    assert min(f_obt) >= 0
    s = sum(f_obt) 
    q = s / sum(f_exp)
    f_exp = [q * x for x in f_exp]
    chisq  = sum ((ni - xi)**2 / xi for ni, xi in zip(f_obt, f_exp))
    p = chisquare_p(chisq, n - 1)
    return chisq, p



######################################################################
# Tests
######################################################################


def do_test_gamma_p():
    """Test functions gamma_p, gamma_q against scipy.special.gammainc"""
    from scipy.special import gammainc
    def gamma_p_testdata():
        for a2 in range(1,30):
            a = 0.5 * a2
            for x10 in range(100):
                 yield a,  0.1 * x10
        yield 100, 300
        yield 100, 5000
    print("Testing incomplete Gamma function")
    for a, x in gamma_p_testdata():
        ref = gammainc(a, x)
        y = gamma_p(a, x)
        assert abs(ref-y) < 1.0e-6, (a, x, ref, y)
        assert abs(y + gamma_q(a,x) -1) < 1.0e-8
    print("passed")


def chisquare_test_distributions():
    """Yield test distributions (f_obt, f_exp) for function chisquare"""
    from random import uniform
    dist_data = [ (3,4,5), (10, 10, 17)]
    
    for n, xmin, xmax in dist_data:
        for i in range(10):
            f_obt = [uniform(xmin, xmax) for i in range(n)]
            f_exp = [uniform(xmin, xmax) for i in range(n)]
            q = sum(f_obt) / sum(f_exp)
            f_exp = [q * x for x in f_exp]
            yield f_obt, f_exp

def do_test_chisquare():
    """Test function chisquare against scipy.stats.chisquare"""
    print("Testing chisquare function")
    from scipy.stats import chisquare as chisquare_ref
    for f_obt, f_exp in chisquare_test_distributions():
         chisq, p = chisquare(f_obt, f_exp)
         chisq_ref, p_ref = chisquare_ref(f_obt, f_exp)
         assert abs(chisq - chisq_ref) < 1.0e-8
         assert abs(p - p_ref)  < 1.0e-6
    print("passed")
          



if __name__ == "__main__":
    do_test_gamma_p()
    do_test_chisquare()

