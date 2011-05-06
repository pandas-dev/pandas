
sparse: pandas/lib/src/sparse.pyx
	-python build_cython.py build_ext --inplace

test: sparse
	-python pandas/lib/tests/test_libsparse.py