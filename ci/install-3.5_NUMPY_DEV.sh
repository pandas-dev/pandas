#!/bin/bash

source activate pandas

echo "install numpy master wheel"

# remove the system installed numpy
pip uninstall numpy -y

# we need these for numpy

# these wheels don't play nice with the conda libgfortran / openblas
# time conda install -n pandas libgfortran openblas || exit 1

# install numpy wheel from master
pip install --pre --upgrade --no-index --timeout=60 --trusted-host travis-dev-wheels.scipy.org -f http://travis-dev-wheels.scipy.org/ numpy scipy

true
