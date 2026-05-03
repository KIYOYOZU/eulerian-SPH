@echo off
setlocal EnableExtensions EnableDelayedExpansion

set "ROOT_DIR=%~dp0..\..\..\..\.."
for %%i in ("%ROOT_DIR%") do set "ROOT_DIR=%%~fi"
set "BUILD_DIR=%ROOT_DIR%\build"
set "TARGET=test_2d_eulerian_flow_around_cylinder_LG"
set "REMAP_TARGET=%TARGET%_remap"
set "CONFIG=Release"
set "VCPKG_ROOT=D:\AAA_postgraduate\SPH\code\vcpkg"
set "CASE_DIR=%~dp0."
for %%i in ("%CASE_DIR%") do set "CASE_DIR=%%~fi"
set "CASE_BIN_DIR=%CASE_DIR%\bin"

cmake --fresh -S "%ROOT_DIR%" -B "%BUILD_DIR%" -A x64 ^
  -DCMAKE_TOOLCHAIN_FILE="%VCPKG_ROOT%\scripts\buildsystems\vcpkg.cmake" ^
  -DSPHINXSYS_2D=ON ^
  -DSPHINXSYS_3D=OFF ^
  -DSPHINXSYS_BUILD_TESTS=ON ^
  -DSPHINXSYS_BUILD_2D_EXAMPLES=ON ^
  -DSPHINXSYS_2D_EXAMPLE_ONLY=%TARGET% ^
  -DSPHINXSYS_BUILD_3D_EXAMPLES=OFF ^
  -DSPHINXSYS_BUILD_PYTHON_INTERFACE=OFF ^
  -DSPHINXSYS_BUILD_UNIT_TESTS=OFF ^
  -DSPHINXSYS_BUILD_EXTRA_SOURCE_AND_TESTS=OFF ^
  -DSPHINXSYS_BUILD_OPTIMIZATION_EXAMPLES=OFF ^
  -DSPHINXSYS_BUILD_MODULES=OFF ^
  -DSPHINXSYS_BUILD_SYCL_TESTS=OFF

if errorlevel 1 exit /b %errorlevel%

cmake --build "%BUILD_DIR%" --config %CONFIG% --target %TARGET% %REMAP_TARGET%
if errorlevel 1 exit /b %errorlevel%

if not exist "%CASE_BIN_DIR%\" mkdir "%CASE_BIN_DIR%"
copy /Y "%BUILD_DIR%\tests\2d_examples\%TARGET%\bin\%CONFIG%\%TARGET%.exe" "%CASE_BIN_DIR%\%TARGET%.exe" >nul
if errorlevel 1 exit /b %errorlevel%
copy /Y "%BUILD_DIR%\tests\2d_examples\%TARGET%\bin\%CONFIG%\%REMAP_TARGET%.exe" "%CASE_BIN_DIR%\%REMAP_TARGET%.exe" >nul
if errorlevel 1 exit /b %errorlevel%

pushd "%CASE_DIR%"
"%CASE_BIN_DIR%\%TARGET%.exe" %*
set "RUN_EXIT_CODE=!ERRORLEVEL!"
popd
exit /b !RUN_EXIT_CODE!
