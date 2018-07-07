tseries: pandas/_libs/lib.pyx pandas/_libs/tslib.pyx pandas/_libs/hashtable.pyx
	python setup.py build_ext --inplace

.PHONY : develop build clean clean_pyc tseries doc

clean:
	-python setup.py clean

clean_pyc:
	-find . -name '*.py[co]' -exec rm {} \;

build: clean_pyc
	python setup.py build_ext --inplace

lint-diff:
	git diff master --name-only -- "*.py" | grep "pandas" | xargs flake8

develop: build
	-python setup.py develop

doc:
	-rm -rf doc/build doc/source/generated
	cd doc; \
	python make.py clean; \
	python make.py html
	python make.py spellcheck
