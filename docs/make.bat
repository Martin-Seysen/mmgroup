REM @ECHO OFF
SETLOCAL

pushd "%~dp0"

REM Command file for Sphinx documentation

SET SPHINXBUILD=python -m sphinx


SET SOURCEDIR=source
SET BUILDDIR=build
SET SPHINXOPTS=-E -a

IF "%~1"=="" GOTO help

REM Check if sphinx-build exists
%SPHINXBUILD% >NUL 2>NUL
IF ERRORLEVEL 9009 (
    ECHO.
    ECHO The 'sphinx-build' command was not found. Make sure you have Sphinx installed,
    ECHO then set the SPHINXBUILD environment variable to point to the full path
    ECHO of the 'sphinx-build' executable. Alternatively, add the Sphinx directory to PATH.
    ECHO.
    ECHO If you don't have Sphinx installed, get it from:
    ECHO https://www.sphinx-doc.org/
    EXIT /B 1
)

%SPHINXBUILD% -M "%~1" "%SOURCEDIR%" "%BUILDDIR%" %SPHINXOPTS% %O%
GOTO end

:help
%SPHINXBUILD% -M help "%SOURCEDIR%" "%BUILDDIR%" %SPHINXOPTS% %O%

:end
popd
ENDLOCAL
