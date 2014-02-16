@echo off
echo "starting 26-64"

setlocal EnableDelayedExpansion
set MSSdk=1
CALL "C:\Program Files\Microsoft SDKs\Windows\v7.0\Bin\SetEnv.cmd" /x64 /release
set DISTUTILS_USE_SDK=1

title 26-64 build
echo "building"
cd "c:\users\Jeff Reback\documents\github\pandas"
C:\python26-64\python.exe setup.py build > build.26-64.log 2>&1

echo "installing"
C:\python26-64\python.exe setup.py bdist --formats=wininst > install.26-64.log 2>&1

echo "testing"
C:\python26-64\scripts\nosetests -A "not slow" build\lib.win-amd64-2.6\pandas > test.26-64.log 2>&1

echo "versions"
cd build\lib.win-amd64-2.6
C:\python26-64\python.exe ../../ci/print_versions.py > ../../versions.26-64.log 2>&1


exit
