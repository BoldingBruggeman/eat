mkdir build
REM mkdir install
cd build
set MKLROOT=%PREFIX%\Library
cmake -DPython3_EXECUTABLE="%PYTHON%" -DCMAKE_BUILD_TYPE=Release %RECIPE_DIR%\.. -DCMAKE_INSTALL_PREFIX=%PREFIX% -DMPI_msmpifec_LIBRARY="%MSMPI_LIB64_CONDA_BACKUP%/msmpifec.lib"
if errorlevel 1 exit 1
cmake --build . --config Release --parallel %CPU_COUNT% --target install
if errorlevel 1 exit 1
cd ..
rmdir /S /Q build
REM mkdir %SCRIPTS%
REM copy install\bin\eat_model_gotm.exe %SCRIPTS%\

