#!/bin/bash

source activate pandas

echo "install numpy master wheel"

# remove the system installed numpy
pip uninstall numpy -y

# install numpy wheel from master
pip install --pre --upgrade --no-index --timeout=60 --trusted-host travis-dev-wheels.scipy.org -f http://travis-dev-wheels.scipy.org/ numpy scipy

true
