#!/bin/bash

echo "inside $0"
# Install Dependencies

# workaround for travis ignoring system_site_packages in travis.yml
rm -f $VIRTUAL_ENV/lib/python$TRAVIS_PYTHON_VERSION/no-global-site-packages.txt

function pip_install {
    NAME="$1"
    if [ x"$2" != x"$2" ]; then
        NAME=$2;
    fi
    echo "travis_fold:begin:$NAME"
    pip install $PIP_ARGS  "$1"
    echo "travis_fold:end:$NAME"
}
function apt_install {
    NAME=$1
    if [ x"$2" != x"!2" ]; then
        NAME=$2;
    fi
    echo "travis_fold:begin:$NAME"
    sudo apt-get install $APT_ARGS "$1"
    echo "travis_fold:end:$NAME"
}

# Hard Deps
pip_install cython
pip_install nose
pip_install python-dateutil
pip_install pytz

# try and get numpy as a binary deb

# numpy is preinstalled on 2.x
# if [ ${TRAVIS_PYTHON_VERSION} == "2.7" ]; then
#     sudo apt-get $APT_ARGS install python-numpy;
# fi

if [ ${TRAVIS_PYTHON_VERSION} == "3.2" ]; then
    apt_install  python3-numpy;
elif [ ${TRAVIS_PYTHON_VERSION} == "3.3" ] || [  x"$LOCALE_OVERRIDE" != x""  ]; then # should be >=3,3
    pip_install numpy==1.7.0
else
     pip_install numpy==1.6.1
fi

# Optional Deps
if [ x"$FULL_DEPS" == x"true" ]; then
    echo "Installing FULL_DEPS"
    if [ ${TRAVIS_PYTHON_VERSION} == "2.7" ]; then
        apt_install python-scipy;
    fi

    if [ ${TRAVIS_PYTHON_VERSION} == "3.2" ]; then
        apt_install $APT_ARGS install python3-scipy;
    fi

    if [ ${TRAVIS_PYTHON_VERSION:0:1} == "2" ]; then
        sudo apt-get $APT_ARGS install libhdf5-serial-dev;
        pip_install numexpr
        pip_install tables
    fi

    pip_install matplotlib;
    pip_install openpyxl;
    pip_install xlwt;
    pip_install xlrd ;
    pip_install 'http://downloads.sourceforge.net/project/pytseries/scikits.timeseries/0.91.3/scikits.timeseries-0.91.3.tar.gz?r=' scikit.timeseries
fi

# if [ x"$VBENCH" == x"true" ]; then
#     if [ ${TRAVIS_PYTHON_VERSION:0:1} == "2" ]; then
#         sudo apt-get $APT_ARGS install libhdf5-serial-dev;
#         pip install numexpr
#         pip install tables
#     fi
#     pip $PIP_ARGS install sqlalchemy git+git://github.com/pydata/vbench.git;
# fi

#build and install pandas
echo "travis_fold:begin:pandas-install"
python setup.py build_ext install
echo "travis_fold:end:pandas-install"

#HACK: pandas is a statsmodels dependency
# so we need to install it after pandas
if [ x"$FULL_DEPS" == x"true" ]; then
    pip_install patsy
    # pick recent 0.5dev dec/2012
    pip_install git+git://github.com/statsmodels/statsmodels@c9062e43b8a5f7385537ca95#egg=statsmodels statsmodels
fi;

# make sure the desired locale is generated
if [ x"$LOCALE_OVERRIDE" != x"" ]; then
    # piggyback this build for plotting tests. oh boy.
    pip_install  matplotlib;

    sudo locale-gen "$LOCALE_OVERRIDE"
fi
