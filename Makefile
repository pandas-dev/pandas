.PHONY : clean develop build clean clean_pyc tseries doc

clean: clean_pyc
	-rm -rf build dist
	-find . -name '*.so' -exec rm -f {} \;

clean_pyc:
	-find . -name '*.pyc' -exec rm -f {} \;

tseries: pandas/lib.pyx pandas/tslib.pyx pandas/hashtable.pyx
	python setup.py build_ext --inplace

sparse: pandas/src/sparse.pyx
	python setup.py build_ext --inplace

build: clean_pyc
	python setup.py build_ext --inplace

develop: build
	-python setup.py develop

doc:
	-rm -rf doc/build
	-rm -rf doc/source/generated
	cd doc; \
	python make.py clean; \
	python make.py html
