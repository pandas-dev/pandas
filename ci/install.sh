#!/bin/bash

# There are 2 distinct pieces that get zipped and cached
# - The venv site-packages dir including the installed dependencies
# - The pandas build artifacts, using the build cache support via
#   scripts/use_build_cache.py
#
# if the user opted in to use the cache and we're on a whitelisted fork
# - if the server doesn't hold a cached version of venv/pandas build,
#   do things the slow way, and put the results on the cache server
#   for the next time.
# -  if the cache files are available, instal some necessaries via apt
#    (no compiling needed), then directly goto script and collect 200$.
#

echo "inside $0"

# Install Dependencie
# as of pip 1.4rc2, wheel files are still being broken regularly, this is a known good
# commit. should revert to pypi when a final release is out
pip install -I git+https://github.com/pypa/pip@42102e9deaea99db08b681d06906c2945f6f95e2#egg=pip
pip Install -I https://bitbucket.org/pypa/setuptools/downloads/setuptools-0.8b6.tar.gz
pip install wheel

# comment this line to disable the fetching of wheel files
PIP_ARGS+=" -I --use-wheel --find-links=https://cache27-pypandas.rhcloud.com/"

# Force virtualenv to accpet system_site_packages
rm -f $VIRTUAL_ENV/lib/python$TRAVIS_PYTHON_VERSION/no-global-site-packages.txt

if [ x"$LOCALE_OVERRIDE" != x"" ]; then
    # make sure the locale is available
    # probably useless, since you would need to relogin
    sudo locale-gen "$LOCALE_OVERRIDE"
fi;

#scipy is not included in the cached venv
if [ x"$FULL_DEPS" == x"true" ] ; then
   # for pytables gets the lib as well
  time sudo apt-get $APT_ARGS install libhdf5-serial-dev

   if [ ${TRAVIS_PYTHON_VERSION:0:1} == "3" ]; then
     time sudo apt-get $APT_ARGS install python3-bs4
   elif [ ${TRAVIS_PYTHON_VERSION:0:1} == "2" ]; then
     time sudo apt-get $APT_ARGS install python-bs4
   fi

   if [ ${TRAVIS_PYTHON_VERSION} == "3.2" ]; then
       time sudo apt-get $APT_ARGS install python3-scipy
   elif [ ${TRAVIS_PYTHON_VERSION} == "2.7" ]; then
       time sudo apt-get $APT_ARGS install python-scipy
   fi
fi

# Hard Deps
time pip install $PIP_ARGS nose python-dateutil pytz>=2013a
time pip install $PIP_ARGS cython==0.19.1

if [ ${TRAVIS_PYTHON_VERSION} == "3.3" ]; then # should be >=3,3
    time pip install $PIP_ARGS numpy==1.7.1
elif [ ${TRAVIS_PYTHON_VERSION} == "3.2" ]; then
    # sudo apt-get $APT_ARGS install python3-numpy; # 1.6.2 or precise
    time pip install $PIP_ARGS numpy==1.6.1
else
    time pip install $PIP_ARGS numpy==1.6.1
fi

# Optional Deps
if [ x"$FULL_DEPS" == x"true" ]; then
    echo "Installing FULL_DEPS"

    if [ ${TRAVIS_PYTHON_VERSION:0:1} == "2" ]; then
        time pip install $PIP_ARGS xlwt
        time pip install $PIP_ARGS bottleneck==0.6.0
        time pip install $PIP_ARGS numexpr==2.1
        time pip install $PIP_ARGS tables==2.3.1
    else
        time pip install $PIP_ARGS numexpr==2.1
        time pip install $PIP_ARGS tables==3.0.0
    fi

    time pip install $PIP_ARGS matplotlib==1.2.1
    time pip install $PIP_ARGS openpyxl
    time pip install $PIP_ARGS xlrd>=0.9.0
    time pip install $PIP_ARGS 'http://downloads.sourceforge.net/project/pytseries/scikits.timeseries/0.91.3/scikits.timeseries-0.91.3.tar.gz?r='
    time pip install $PIP_ARGS patsy
    time pip install $PIP_ARGS html5lib

    if [ ${TRAVIS_PYTHON_VERSION:0:1} == "3" ]; then
      time sudo apt-get $APT_ARGS remove python3-lxml
    elif [ ${TRAVIS_PYTHON_VERSION:0:1} == "2" ]; then
      time sudo apt-get $APT_ARGS remove python-lxml
    fi

    # fool statsmodels into thinking pandas was already installed
    # so it won't refuse to install itself.

    mkdir  $SITE_PKG_DIR/pandas
    touch $SITE_PKG_DIR/pandas/__init__.py
    echo "version='0.10.0-phony'" >  $SITE_PKG_DIR/pandas/version.py
    time pip install $PIP_ARGS git+git://github.com/statsmodels/statsmodels@c9062e43b8a5f7385537ca95#egg=statsmodels

    rm -Rf $SITE_PKG_DIR/pandas # scrub phoney pandas
fi

# build pandas
time python setup.py build_ext install

true
