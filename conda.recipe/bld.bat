mkdir build
REM mkdir install
cd build

set MKLROOT=%PREFIX%\Library

set EXTERNAL_MODEL_DIR=%RECIPE_DIR%\..\extern
set EXTERNAL_INSTITUTES=ogs
set MODEL_DIRS=-DFABM_ERSEM_BASE=%EXTERNAL_MODEL_DIR%\ersem -DFABM_OGS_BASE=%EXTERNAL_MODEL_DIR%\ogs -DFABM_PISCES_BASE=%EXTERNAL_MODEL_DIR%\pisces

cmake -DFABM_EXTRA_INSTITUTES=%EXTERNAL_INSTITUTES% %MODEL_DIRS% -DPython3_EXECUTABLE="%PYTHON%" -DCMAKE_BUILD_TYPE=Release %RECIPE_DIR%\.. -DCMAKE_INSTALL_PREFIX=%PREFIX% -DMPI_msmpifec_LIBRARY="%MSMPI_LIB64_CONDA_BACKUP%/msmpifec.lib"
if errorlevel 1 exit 1
cmake --build . --config Release --parallel %CPU_COUNT% --target install
if errorlevel 1 exit 1
cd ..
rmdir /S /Q build
