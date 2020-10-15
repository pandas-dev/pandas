git bisect start ${1} ${2}
git bisect run bash -c "python setup.py build_ext --inplace -j 4; /
    python -m pip install -e . --no-build-isolation --no-use-pep517; /
    python ${3}"
