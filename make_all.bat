rem *******************************************************
rem Create external modules and documentation
rem *******************************************************

python setup.py   build_ext --nprocesses 16 --inplace -f -cmsvc
if errorlevel 1 goto :abort

@echo off
pause Press any key to exit
goto :done 

:abort:
@echo off
pause Press any key to abort
:done
