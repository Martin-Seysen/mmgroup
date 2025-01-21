"""Register custom markers for pytest.

This module registers the custom markers for pytest so 
that we may invoke:

pytest --pyargs mmgroup [options]

without 'PytestUnknownMarkWarning' warnings. 
This tests the **installed** mmgroup package.

This registering is usually done in the configuration 
file pytest.ini in the root directory, see:

https://pytest.org/en/7.4.x/reference/customize.html#pytest-ini

Alternatively, we may register these markers in file conftest.py.
Here file conftest.py may be loacated in the  installed
python package mmgroup, see:

https://pytest.org/en/7.4.x/how-to/writing_plugins.html#registering-custom-markers

"""


# Registration of markers in the style used by file pytest.ini 
markers = r"""
   auto_group: test groups derived from mmgroup.structures.auto_group
   axes:       test modules dealing with operations on 2A axes
   bench:      benchmark
   bimm:       test modules implementing the BiMonster
   bitfunc:    test for the bitfunctions module
   build:      tests to be done after build process (very few tests)
   compiler:   test requires a C compiler 
   demo:       test demonstration code for reduction algorithm
   extremely_slow:  marks tests as even slower than very slow
   hadamard:   test exection of code generated for hadamard matrices
   involution: test module involution.c
   gen_xi:     test for functions gen_XXX in the generators module
   general:    test general group operations and union-find algorithm
   mat24:      test for the mat24 module
   mm:         test high-level group operation
   mm_op:      test group operation on vector space in mmgroup.mm_space
   mm_op_crt:  test group operation on vector space in mmgroup.mm_crt_space
   mmgroup:    test group operations in mmgroup.mm_group
   orders:     test function for orders and characters ins monster group
   ploop:      test classes for Golay code, Parker loop, AutPLoop.
   qstate:     test (subclasses of) class  mmgroup.clifford12.QState12 
   slow:       marks tests as slow (deselect with '-m "not slow"')
   space:      test spaces derived from mmgroup.structures.abstract_mm_rep_space
   user:       interaction tests that should be done by the user
   very_slow:  marks tests as very slow 
   xsp2co1:    test modules xsp2co1.c and xsp2co1_elem.c
"""

# Hook for extending pytest configuration
def pytest_configure(config):
    for line in markers.split("\n"):
        if len(line) and not line.isspace():
            config.addinivalue_line("markers", line.strip())
