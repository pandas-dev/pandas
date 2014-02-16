@echo off
echo "starting 26-32"


title 26-32 build
echo "building"
cd "c:\users\Jeff Reback\documents\github\pandas"
C:\python26-32\python.exe setup.py build > build.26-32.log 2>&1

echo "installing"
C:\python26-32\python.exe setup.py bdist --formats=wininst > install.26-32.log 2>&1

echo "testing"
C:\python26-32\scripts\nosetests -A "not slow" build\lib.win32-2.6\pandas > test.26-32.log 2>&1

echo "versions"
cd build\lib.win32-2.6
C:\python26-32\python.exe ../../ci/print_versions.py > ../../versions.26-32.log 2>&1


exit
