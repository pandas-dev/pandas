@echo off
echo "starting 27-32"

setlocal EnableDelayedExpansion
set MSSdk=1
CALL "C:\Program Files\Microsoft SDKs\Windows\v7.0\Bin\SetEnv.cmd" /x86 /release
set DISTUTILS_USE_SDK=1

title 27-32 build
echo "building"
cd "c:\users\Jeff Reback\documents\github\pandas"
C:\python27-32\python.exe setup.py build > build.27-32.log 2>&1

title "installing"
C:\python27-32\python.exe setup.py bdist --formats=wininst > install.27-32.log 2>&1

echo "testing"
C:\python27-32\scripts\nosetests -A "not slow" build\lib.win32-2.7\pandas > test.27-32.log 2>&1

echo "versions"
cd build\lib.win32-2.7
C:\python27-32\python.exe ../../ci/print_versions.py > ../../versions.27-32.log 2>&1

exit

