project('mmgroup', 'cython', 'c')

buildtype = 'release'

prog_python = find_program('python3')
py = import('python').find_installation()
dep_py = py.dependency()


C_DIR = meson.build_root() + '/c_files'
incdir = include_directories('c_files')

MOCKUP = ''

GENERATE_START = '''
 -v
 @0@
 --py-path @SOURCE_ROOT@/src
           @OUTDIR@
 --out-dir @OUTDIR@/c_files
 --out-pxd-dir @OUTDIR@/c_files
'''.format(MOCKUP, C_DIR)


GEN_CODE =  [prog_python, '@SOURCE_ROOT@/generate_code.py'] 


MAT24_SOURCES = '''
   mat24_functions.c
   mat24_random.c
'''

MAT24_GENERATE = GENERATE_START + '''
 --dll MAT24
 --source-path @SOURCE_ROOT@/src/mmgroup/dev/mat24
 --tables mmgroup.dev.mat24.mat24_ref 
 --sources mat24_functions.h
 --sources
''' + MAT24_SOURCES + '''
 --pxd  mat24_functions.pxd
 --pyx  mat24fast.pyx
'''



MAT24_SOURCE_PATHS = []
foreach f : MAT24_SOURCES.split()
    MAT24_SOURCE_PATHS += [ C_DIR + '/' + f]
endforeach

message('MAT24_SOURCE_PATHS =', MAT24_SOURCE_PATHS)


mat24_c_target = custom_target('target_mat24',
  output : ['mat24fast.pyx'] + MAT24_SOURCES.split(),
  input : [],
  command : GEN_CODE + MAT24_GENERATE.split(),
  install : true, 
  install_dir : C_DIR, 
)

#mat24_dep = declare_dependency (sources: mat24_c_target)

message( 'mat24_c_target path',  mat24_c_target[0].full_path())



py.extension_module(
    'mat24',
     sources : mat24_c_target,
     dependencies : [dep_py],
     include_directories : incdir,
     install_dir: C_DIR + '/mmgroup',
     native : true, 
)