mkdir build
REM mkdir install
cd build
cmake -DPython3_EXECUTABLE="%PYTHON%" -DCMAKE_BUILD_TYPE=Release %RECIPE_DIR%\.. -DFABM_BASE=%RECIPE_DIR%/../extern/fabm -DCMAKE_INSTALL_PREFIX=%PREFIX%
if errorlevel 1 exit 1
cmake --build . --config Release --parallel %CPU_COUNT% --target install
if errorlevel 1 exit 1
cd ..
rmdir /S /Q build
REM mkdir %SCRIPTS%
REM copy install\bin\eat_model_gotm.exe %SCRIPTS%\

