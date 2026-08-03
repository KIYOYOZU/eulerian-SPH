@echo off
setlocal

set "ROOT_DIR=%~dp0..\..\.."
for %%i in ("%ROOT_DIR%") do set "ROOT_DIR=%%~fi"
set "BUILD_DIR=%ROOT_DIR%\build"
set "TARGET=test_2d_eulerian_supersonic_flow_new_BC"
set "CONFIG=Release"
set "VCPKG_ROOT=D:\AAA_postgraduate\SPH\code\vcpkg"
set "CASE_BIN_DIR=%~dp0bin"

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

cmake --build "%BUILD_DIR%" --config %CONFIG% --target %TARGET%
if errorlevel 1 exit /b %errorlevel%

if not exist "%CASE_BIN_DIR%\" mkdir "%CASE_BIN_DIR%"
copy /Y "%BUILD_DIR%\tests\2d_examples\%TARGET%\bin\%CONFIG%\%TARGET%.exe" "%CASE_BIN_DIR%\%TARGET%.exe" >nul
if errorlevel 1 exit /b %errorlevel%

pushd "%~dp0"
"%CASE_BIN_DIR%\%TARGET%.exe" %*
popd
