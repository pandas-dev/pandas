#!/bin/bash

echo "inside $0"
# Install Dependencies

# Hard Deps
pip install $PIP_ARGS --use-mirrors cython nose python-dateutil pytz

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

    if [ ${TRAVIS_PYTHON_VERSION:0:1} == "2" ]; then
        sudo apt-get $APT_ARGS install libhdf5-serial-dev;
        pip install numexpr
        pip install tables
    fi

    pip install $PIP_ARGS --use-mirrors openpyxl matplotlib;
    pip install $PIP_ARGS --use-mirrors xlrd xlwt;
    pip install $PIP_ARGS 'http://downloads.sourceforge.net/project/pytseries/scikits.timeseries/0.91.3/scikits.timeseries-0.91.3.tar.gz?r='
fi

if [ x"$VBENCH" == x"true" ]; then
    if [ ${TRAVIS_PYTHON_VERSION:0:1} == "2" ]; then
        sudo apt-get $APT_ARGS install libhdf5-serial-dev;
        pip install numexpr
        pip install tables
    fi
    pip $PIP_ARGS install sqlalchemy git+git://github.com/pydata/vbench.git;
fi

#build and install pandas
python setup.py build_ext install

#HACK: pandas is a statsmodels dependency
# so we need to install it after pandas
if [ x"$FULL_DEPS" == x"true" ]; then
    pip install patsy
    # pick recent 0.5dev dec/2012
    pip install git+git://github.com/statsmodels/statsmodels@c9062e43b8a5f7385537ca95#egg=statsmodels
fi;
