tseries: pandas/lib.pyx pandas/tslib.pyx pandas/hashtable.pyx
	python setup.py build_ext --inplace

.PHONY : develop build clean clean_pyc tseries doc

clean:
	-python setup.py clean

clean_pyc:
	-find . -name '*.py[co]' -exec rm {} \;

sparse: pandas/src/sparse.pyx
	python setup.py build_ext --inplace

build: clean_pyc
	python setup.py build_ext --inplace

develop: build
	-python setup.py develop

doc:
	-rm -rf doc/build doc/source/generated
	cd doc; \
	python make.py clean; \
	python make.py html
