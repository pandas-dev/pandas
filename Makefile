.PHONY : develop build clean clean_pyc doc lint-diff black

all: develop

clean:
	-python setup.py clean

clean_pyc:
	-find . -name '*.py[co]' -exec rm {} \;

build: clean_pyc
	python setup.py build_ext --inplace

lint-diff:
	git diff upstream/master --name-only -- "*.py" | xargs flake8

black:
	black .

develop: build
	python -m pip install --no-build-isolation -e .

doc:
	-rm -rf doc/build doc/source/generated
	cd doc; \
	python make.py clean; \
	python make.py html
