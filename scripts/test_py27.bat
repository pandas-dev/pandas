SET PATH=C:\MinGW\bin;C:\Python27;C:\Python27\Scripts;%PATH%

python setup.py clean
python setup.py build_ext -c mingw32 --inplace

nosetests pandas