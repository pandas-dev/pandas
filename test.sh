#!/bin/sh
coverage erase
nosetests -w pandas/core --with-coverage --cover-package=pandas.core --pdb-failure --pdb
# nosetests -w pandas/stats --with-coverage --cover-package=pandas.stats
# coverage run runtests.py