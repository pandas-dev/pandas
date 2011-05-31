clean:
	-rm -rf build dist

tseries: pandas/lib/src/tseries.pyx
	touch pandas/lib/src/tseries.pyx
	python build_cython.py build_ext --inplace

sparse: pandas/lib/src/sparse.pyx
	-python build_cython.py build_ext --inplace

test: sparse
	-python pandas/lib/tests/test_libsparse.py