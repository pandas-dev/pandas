SET PATH=C:\MinGW\bin;C:\Python27;C:\Python27\Scripts;%PATH%
del pandas\lib\tseries.pyd
python setup.py build_ext -c mingw32 --inplace
nosetests pandas