#!/bin/sh
#nosetests -w pandas/core --with-coverage --cover-package=pandas.core --pdb-failure
coverage erase
nosetests -w pandas/core --with-coverage --cover-package=pandas.core
# coverage run runtests.py