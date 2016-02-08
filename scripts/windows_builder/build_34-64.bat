@echo off
echo "starting 34-64"

setlocal EnableDelayedExpansion
set MSSdk=1
CALL "C:\Program Files\Microsoft SDKs\Windows\v7.1\Bin\SetEnv.cmd" /x64 /release
set DISTUTILS_USE_SDK=1

title 34-64 build
echo "building"
cd "c:\users\Jeff Reback\documents\github\pandas"
C:\python34-64\python.exe setup.py build > build.34-64.log 2>&1

echo "installing"
C:\python34-64\python.exe setup.py bdist --formats=wininst > install.34-64.log 2>&1

echo "testing"
C:\python34-64\scripts\nosetests -A "not slow" build\lib.win-amd64-3.4\pandas > test.34-64.log 2>&1

echo "versions"
cd build\lib.win-amd64-3.4
C:\python34-64\python.exe ../../ci/print_versions.py > ../../versions.34-64.log 2>&1

exit



