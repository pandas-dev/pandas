#!/bin/sh
coverage erase
# nosetests -w pandas --with-coverage --cover-package=pandas --pdb-failure --pdb
nosetests -w pandas/core --with-coverage --cover-package=pandas.core --pdb-failure --pdb
# nosetests -w pandas/stats --with-coverage --cover-package=pandas.stats
# coverage run runtests.py