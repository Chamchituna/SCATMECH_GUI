@echo off
setlocal EnableExtensions DisableDelayedExpansion
cd /d "%~dp0"
set "PYTHON_EXE="
set "PYLAUNCHER_DRYRUN="
set "PYLAUNCHER_ALLOW_INSTALL="
set "PY_CHECK=import sys,struct,platform; sys.exit(0 if sys.version_info[:2] == (3,12) and struct.calcsize('P') == 8 and platform.machine().lower() in ('amd64','x86_64') and platform.python_implementation() == 'CPython' else 1)"

if not "%~1"=="" goto explicit_python

py.exe -3.12 -c "%PY_CHECK%" >nul 2>nul
if "%ERRORLEVEL%"=="0" goto run_launcher

call :try_python "%LOCALAPPDATA%\Programs\Python\Python312\python.exe"
if defined PYTHON_EXE goto run_python
call :try_python "%ProgramFiles%\Python312\python.exe"
if defined PYTHON_EXE goto run_python
call :try_python "%SystemDrive%\Python312\python.exe"
if defined PYTHON_EXE goto run_python

for /f "delims=" %%P in ('where python.exe 2^>nul') do call :try_python "%%P"
if defined PYTHON_EXE goto run_python
for /f "delims=" %%P in ('where python3.exe 2^>nul') do call :try_python "%%P"
if defined PYTHON_EXE goto run_python
goto missing_python

:explicit_python
call :try_python "%~1"
if defined PYTHON_EXE goto run_python
echo The supplied executable is not a working CPython 3.12 x64 installation.
goto missing_python

:try_python
if defined PYTHON_EXE exit /b 0
if not exist "%~1" exit /b 1
"%~1" -c "%PY_CHECK%" >nul 2>nul
if not "%ERRORLEVEL%"=="0" exit /b 1
set "PYTHON_EXE=%~1"
exit /b 0

:run_launcher
echo Using Python 3.12 x64 through the Python launcher.
py.exe -3.12 build_windows.py
goto finished

:run_python
echo Using "%PYTHON_EXE%"
"%PYTHON_EXE%" build_windows.py
goto finished

:missing_python
echo.
echo Python 3.12 x64 could not be found. The build has not started.
echo The Microsoft Store shortcut is not a Python installation.
echo.
echo 1. Download and run the official 64-bit installer:
echo    https://www.python.org/ftp/python/3.12.10/python-3.12.10-amd64.exe
echo 2. Select "Add python.exe to PATH" and keep pip and the launcher enabled.
echo 3. Close this window, then double-click BUILD_WINDOWS.cmd again.
echo.
echo If Python is installed in a custom location, run:
echo    BUILD_WINDOWS.cmd "C:\path\to\Python312\python.exe"
echo See README.md for details. No build log is created at this stage.
set "BUILD_RESULT=2"
goto end

:finished
set "BUILD_RESULT=%ERRORLEVEL%"
echo.
if not "%BUILD_RESULT%"=="0" echo Build failed. See the message above; if created, the log is build\build.log.
:end
pause
exit /b %BUILD_RESULT%
