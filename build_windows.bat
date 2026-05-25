@echo off
setlocal
cd /d "%~dp0"

echo === ChromaView Windows build ===
echo.

REM Regenerate icon (requires Pillow, already in venv)
echo [1/3] Generating icon...
call .venv\Scripts\python.exe make_icon.py
if errorlevel 1 ( echo ERROR: icon generation failed & pause & exit /b 1 )

REM Run PyInstaller
echo.
echo [2/3] Running PyInstaller (this takes a minute)...
call .venv\Scripts\pyinstaller.exe --noconfirm ChromaView.spec
if errorlevel 1 ( echo ERROR: PyInstaller failed & pause & exit /b 1 )

echo.
echo [3/3] Done!
echo.
echo Executable: %~dp0dist\ChromaView\ChromaView.exe
echo.
pause
