clean:
	-rm -rf build dist

tseries: pandas/src/tseries.pyx
	python setup.py build_ext --inplace

sparse: pandas/src/sparse.pyx
	-python setup.py build_ext --inplace

test: sparse
	-python pandas/tests/test_libsparse.py