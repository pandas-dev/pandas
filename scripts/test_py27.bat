SET PATH=C:\MinGW\bin;C:\Python27;C:\Python27\Scripts;%PATH%
del pandas\_tseries.pyd
del pandas\_sparse.pyd
python setup.py build_ext -c mingw32 --inplace
nosetests pandas