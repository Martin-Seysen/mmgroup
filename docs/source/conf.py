# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import subprocess
sys.path.insert(0, os.path.abspath(os.path.join('..', '..')))
sys.path.insert(0, os.path.abspath(os.path.join('..', '..', 'src')))


# -- Check if we are in readthedocs -------------------------------------------

on_readthedocs = os.environ.get('READTHEDOCS') == 'True'


# -- Creating C filesin in readthedocs ----------------------------------------

SETUP_DIR = os.path.abspath(os.path.join('..', '..'))
C_DIR = os.path.join(SETUP_DIR, "src", "mmgroup", "dev", "c_files")
DOXYGEN_DIR = os.path.abspath(os.path.join('..', 'doxygen'))




# -- Project information -----------------------------------------------------

project = 'mmgroup'
copyright = '2020, Martin Seysen'
author = 'Martin Seysen'

# The full version, including alpha/beta/rc tags
release = '0.0.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
   'sphinx.ext.autodoc',
   'sphinxcontrib.bibtex',
    'breathe',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
#html_static_path = ['_static']

# -- Options for latex output ------------------------------------------------
# The following option should remove excessive blank pages.

latex_elements = {
  'extraclassoptions': 'openany,oneside'
}



# -- Mock up Cython extension modules ----------------------------------------

# We will mock up most Cython extensions when generating the documentation
# with Sphinx. This means that these extensions need not be present when
# the documentation is generated. For background, see:
# https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html#confval-autodoc_mock_imports

autodoc_mock_imports = [
  #  "mmgroup.mat24",   # We cannot mock up this extension, see below
    "mmgroup.mat24_xi",
    "mmgroup.clifford12",
    "mmgroup.mm",
    "mmgroup.mm3",
    "mmgroup.mm7",
    "mmgroup.mm15",
    "mmgroup.mm31",
    "mmgroup.mm63",
    "mmgroup.mm127",
    "mmgroup.mm255",
    "blah",
]

# For documentation, we cannot mock up the 'mmgroup.mat24' extension,
# since we need some data from that extension.
# Whenever needed, we will import a pure python substitute for the
# 'mmgroup.mat24' extension if the original extension has not been found.
# Class 'Mat24' in module  'mmgroup.dev.mat24.mat24_ref' is the
# best substitute available for 'mmgroup.mat24'.





# -- Workaround for a certain LaTeX error when building the .pdf file --------

# When using nested lists in Sphinx, the following error occurs when building
# the .pdf file:
# ! LaTeX Error: Too deeply nested.
# Here we implement a workaround for this problem, which has been taken from:
# https://stackoverflow.com/questions/28454217/how-to-avoid-the-too-deeply-nested-error-when-creating-pdfs-with-sphinx


latex_elements = {
# Additional stuff for the LaTeX preamble.
'preamble': r'\usepackage{enumitem}\setlistdepth{99}',
}




# -- Breathe Configuration --------

breathe_projects = {
   'mmgroup' : os.path.join(DOXYGEN_DIR, 'xml')
}

breathe_default_project = "mmgroup"

# This is not nice, but C file nmaes must be known in advance
C_FILES = [
"clifford12.h",
"mat24_functions.c",
"mat24_functions.h",
"mat24_xi_functions.c",
"mat24_xi_functions.h",
"mm127_op_misc.c",
"mm127_op_pi.c",
"mm127_op_t.c",
"mm127_op_word.c",
"mm127_op_xi.c",
"mm127_op_xy.c",
"mm15_op_misc.c",
"mm15_op_pi.c",
"mm15_op_t.c",
"mm15_op_word.c",
"mm15_op_xi.c",
"mm15_op_xy.c",
"mm255_op_misc.c",
"mm255_op_pi.c",
"mm255_op_t.c",
"mm255_op_word.c",
"mm255_op_xi.c",
"mm255_op_xy.c",
"mm31_op_misc.c",
"mm31_op_pi.c",
"mm31_op_t.c",
"mm31_op_word.c",
"mm31_op_xi.c",
"mm31_op_xy.c",
"mm3_op_misc.c",
"mm3_op_pi.c",
"mm3_op_t.c",
"mm3_op_word.c",
"mm3_op_xi.c",
"mm3_op_xy.c",
"mm7_op_misc.c",
"mm7_op_pi.c",
"mm7_op_t.c",
"mm7_op_word.c",
"mm7_op_xi.c",
"mm7_op_xy.c",
"mm_aux.c",
"mm_basics.h",
"mm_crt.c",
"mm_group_n.c",
"mm_group_word.c",
"mm_op127.h",
"mm_op15.h",
"mm_op255.h",
"mm_op3.h",
"mm_op31.h",
"mm_op7.h",
"mm_random.c",
"mm_tables.c",
"mm_tables_xi.c",
"qmatrix12.c",
"qstate12.c",
"xsp2co1.c",
]

breathe_projects_source = {
    "mmgroup" :  (C_DIR, C_FILES) # C_FILES not yet known
}


# -- Call generate_doxygen_xml -----------------------------------------------

def generate_doxygen_xml(app):
    """Run the doxygen make commands if we're on the ReadTheDocs server"""

    print(" \nStarting generate_doxygen_xml ...")
    if on_readthedocs:
        print("\nGenerating C files ...")
        subprocess.check_call([sys.executable, "setup.py", "build_ext"], 
            cwd=SETUP_DIR)
        print("C files have been generated\n")

    print("Doxygen Directory = ", DOXYGEN_DIR)
    print("C Directory =", C_DIR, ":")
    C_FILES_FOUND = [f for f in os.listdir(C_DIR) 
        if os.path.splitext(f)[1] in [".c", ".h"]]
    print(C_FILES_FOUND, "\n")
    subprocess.check_call("doxygen", cwd = DOXYGEN_DIR)
    print("End of generate_doxygen_xml\n")



def setup(app):
    # For readthedocs:
    # Add hook for building doxygen xml when needed
    app.connect("builder-inited", generate_doxygen_xml)
    pass

