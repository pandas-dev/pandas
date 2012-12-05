#!/bin/bash

echo "inside $0"
# Install Dependencies

# Hard Deps
pip install $PIP_ARGS --use-mirrors cython nose python-dateutil

# try and get numpy as a binary deb

# numpy is preinstalled on 2.x
# if [ ${TRAVIS_PYTHON_VERSION} == "2.7" ]; then
#     sudo apt-get $APT_ARGS install python-numpy;
# fi

if [ ${TRAVIS_PYTHON_VERSION} == "3.2" ]; then
    sudo apt-get $APT_ARGS install python3-numpy;
fi

# or else, get it with pip and compile it
if [ ${TRAVIS_PYTHON_VERSION:0:1} == "2" ] || \
   [ ${TRAVIS_PYTHON_VERSION}     == "3.1" ] || \
   [ ${TRAVIS_PYTHON_VERSION}     == "3.2" ]; then
     pip $PIP_ARGS install numpy; #https://github.com/y-p/numpy/archive/1.6.2_with_travis_fix.tar.gz;
else
    pip $PIP_ARGS install https://github.com/numpy/numpy/archive/v1.7.0b2.tar.gz;
fi

# Optional Deps
if [ x"$FULL_DEPS" == x"true" ]; then
    echo "Installing FULL_DEPS"
    if [ ${TRAVIS_PYTHON_VERSION} == "2.7" ]; then
        sudo apt-get $APT_ARGS install python-scipy;
    fi

    if [ ${TRAVIS_PYTHON_VERSION} == "3.2" ]; then
        sudo apt-get $APT_ARGS install python3-scipy;
    fi

    pip install $PIP_ARGS --use-mirrors openpyxl pytz matplotlib;
    pip install $PIP_ARGS --use-mirrors xlrd xlwt;
fi

if [ x"$VBENCH" == x"true" ]; then
    pip $PIP_ARGS install sqlalchemy git+git://github.com/pydata/vbench.git;
fi
