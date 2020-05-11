rem *******************************************************
rem Create documentation
rem *******************************************************

cd docs
call make html
if errorlevel 1 goto :abort_docs
call make latexpdf
if errorlevel 1 goto :abort_docs
cd ..

@echo off
pause Press any key to exit
goto :done 

:abort_docs:
cd ..
:abort:
@echo off
pause Press any key to abort
exit 1
:done
