clean:
	-rm -rf build dist

tseries: pandas/lib.pyx pandas/tslib.pyx pandas/hashtable.pyx
	python setup.py build_ext --inplace

sparse: pandas/src/sparse.pyx
	-python setup.py build_ext --inplace

test: sparse
	-python pandas/tests/test_libsparse.py