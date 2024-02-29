set REPO_DIR=%~dp0
set WORK_DIR=%TEMP%\build-eat-%RANDOM%
cmake -B "%WORK_DIR%" -S "%REPO_DIR%" %*
if errorlevel 1 exit /b 1
cmake --build "%WORK_DIR%" --parallel 4 --target install --config Release
if errorlevel 1 exit /b 1
rmdir /S /Q "%WORK_DIR%"
