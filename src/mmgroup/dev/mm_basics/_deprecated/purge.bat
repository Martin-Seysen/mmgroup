del mm*.c
del mm*.h
del mm*.o
del mm*.pxd
del mm*.pyx
del *.pyc
rd /s /q build
cd codegen_sample
call purge.bat
cd..
cd cython_sample
call purge.bat
cd..





