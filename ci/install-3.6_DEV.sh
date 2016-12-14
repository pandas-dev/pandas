#!/bin/bash

echo "install 3.6 dev"

conda config --set add_pip_as_python_dependency false
conda create -n pandas python=3.6 -c conda-forge/label/prerelease

source activate pandas

# ensure we have pip
python -m ensurepip

# install deps
pip3.6 install nose cython numpy pytz python-dateutil

true
