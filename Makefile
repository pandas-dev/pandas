.PHONY : develop build clean clean_pyc doc lint-diff black check

all: develop

clean:
	-python3 setup.py clean

clean_pyc:
	-find . -name '*.py[co]' -exec rm {} \;

build: clean_pyc
	python3 setup.py build_ext --inplace

lint-diff:
	git diff upstream/master --name-only -- "*.py" | xargs flake8

black:
	black .

develop: build
	python3 -m pip install --no-build-isolation -e .

doc:
	-rm -rf doc/build doc/source/generated
	cd doc; \
	python3 make.py clean; \
	python3 make.py html

check:
	python3 scripts/validate_unwanted_patterns.py \
		--validation-type="private_function_across_module" \
		--included-file-extensions="py" \
		--excluded-file-paths=pandas/tests,asv_bench/,pandas/_vendored \
		pandas/

	python3 scripts/validate_unwanted_patterns.py \
		--validation-type="private_import_across_module" \
		--included-file-extensions="py" \
		--excluded-file-paths=pandas/tests,asv_bench/,pandas/_vendored,doc/
		pandas/
