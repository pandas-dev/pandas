#!/bin/bash

echo "install 3.6 dev"

conda config --set add_pip_as_python_dependency false
conda create -n pandas python=3.6 -c conda-forge/label/prerelease

source activate pandas

# ensure we have pip
python -m ensurepip
pip3.6 install nose

# build cython
git clone https://github.com/cython/cython.git
cd cython
git checkout 0.25.1
python setup.py install
cd ..

# remove the system installed numpy
pip3.6 uninstall numpy -y

# install deps
pip3.6 install numpy
pip3.6 install pytz
pip3.6 install python-dateutil

true
