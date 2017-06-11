#!/bin/bash

source activate pandas

echo "install numpy master wheel"

# remove the system installed numpy
pip uninstall numpy -y

# install numpy wheel from master
PRE_WHEELS="https://7933911d6844c6c53a7d-47bd50c35cd79bd838daf386af554a83.ssl.cf2.rackcdn.com"
pip install --pre --upgrade --timeout=60 -f $PRE_WHEELS numpy scipy

# install dateutil from master
pip install -U git+git://github.com/dateutil/dateutil.git

true
