#!/bin/sh

python setup.py clean
python setup.py build_ext --inplace
coverage erase
# nosetests pandas/tests/test_index.py --with-coverage --cover-package=pandas.core --pdb-failure --pdb
#nosetests -w pandas --with-coverage --cover-package=pandas --pdb-failure --pdb #--cover-inclusive
nosetests -w pandas --with-coverage --cover-package=pandas $* #--cover-inclusive
# nosetests -w pandas/io --with-coverage --cover-package=pandas.io --pdb-failure --pdb
# nosetests -w pandas/core --with-coverage --cover-package=pandas.core --pdb-failure --pdb
# nosetests -w pandas/stats --with-coverage --cover-package=pandas.stats
# coverage run runtests.py
