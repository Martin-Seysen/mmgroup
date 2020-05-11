@echo on
python codegen.py
python setup.py build_ext --inplace -c mingw32
@echo off
if errorlevel 1 goto abort
@echo on
python setup.py clean
del mm_op.c
cd ..\test_mm
call testit.bat %1
if errorlevel 1 goto abort
cd ..\mm

@echo off
if "%1" == "b" goto done
pause "Press any key to exit!"
:done
exit /b 0
:abort
pause "Press any key to abort!"
exit 1

