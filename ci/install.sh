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

function edit_init() {
    if [ -n "$LOCALE_OVERRIDE" ]; then
        pandas_dir=pandas
        echo "Adding locale to the first line of $pandas_dir/__init__.py"
        rm -f $pandas_dir/__init__.pyc
        sedc="1iimport locale; locale.setlocale(locale.LC_ALL, '$LOCALE_OVERRIDE')"
        sed -i "$sedc" $pandas_dir/__init__.py
        echo "First line of $pandas_dir/__init__.py"
        head $pandas_dir/__init__.py
    fi
}

echo "inside $0"

# Install Dependencies
# as of pip 1.4rc2, wheel files are still being broken regularly, this is a known good
# commit. should revert to pypi when a final release is out
pip install -I git+https://github.com/pypa/pip@42102e9deaea99db08b681d06906c2945f6f95e2#egg=pip
pv="${TRAVIS_PYTHON_VERSION:0:1}"
[ "$pv" == "2" ] && pv=""

pip install -I -U setuptools
pip install wheel

# comment this line to disable the fetching of wheel files
PIP_ARGS+=" -I --use-wheel --find-links=http://cache27diy-cpycloud.rhcloud.com/${TRAVIS_PYTHON_VERSION}${JOB_TAG}/"

# Force virtualenv to accpet system_site_packages
rm -f $VIRTUAL_ENV/lib/python$TRAVIS_PYTHON_VERSION/no-global-site-packages.txt


if [ -n "$LOCALE_OVERRIDE" ]; then
    # make sure the locale is available
    # probably useless, since you would need to relogin
    sudo locale-gen "$LOCALE_OVERRIDE"
fi


# show-skipped is working at this particular commit
time pip install git+git://github.com/cpcloud/nose-show-skipped.git@fa4ff84e53c09247753a155b428c1bf2c69cb6c3
time pip install $PIP_ARGS -r ci/requirements-${TRAVIS_PYTHON_VERSION}${JOB_TAG}.txt
time sudo apt-get install libatlas-base-dev gfortran


# install gui for clipboard testing
if [ -n "$GUI" ]; then
    echo "Using GUI clipboard: $GUI"
    [ -n "$pv" ] && py="py"
    time sudo apt-get $APT_ARGS install python${pv}-${py}${GUI}
fi


# install a clipboard
if [ -n "$CLIPBOARD" ]; then
    echo "Using clipboard: $CLIPBOARD"
    time sudo apt-get $APT_ARGS install $CLIPBOARD
fi


# Optional Deps
if [ x"$FULL_DEPS" == x"true" ]; then
    echo "Installing FULL_DEPS"
    # for pytables gets the lib as well
    time sudo apt-get $APT_ARGS install libhdf5-serial-dev
fi


edit_init

# build pandas
time python setup.py build_ext install

true
