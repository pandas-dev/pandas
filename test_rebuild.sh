#!/bin/sh

python setup.py clean
python setup.py build_ext --inplace
coverage erase
pytest pandas --cov=pandas
