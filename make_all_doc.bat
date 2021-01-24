call msvc.bat

python setup.py  build_ext --inplace -f -cmsvc
if errorlevel 1 goto :abort

call make_doc
if errorlevel 1 goto :abort

:abort:
@echo off
pause Press any key to abort
:done
