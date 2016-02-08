#!/bin/bash

# This script is meant to run on a mint precise64 VM.
# The generated wheel files should be compatible
# with travis-ci as of 07/2013.
#
# Runtime can be up to an hour or more.

echo "Building wheels..."

# print a trace for everything; RTFM
set -x

# install and update some basics
apt-get update
apt-get install python-software-properties git -y
apt-add-repository ppa:fkrull/deadsnakes -y
apt-get update

# install some deps and virtualenv
apt-get install python-pip libfreetype6-dev libpng12-dev libhdf5-serial-dev \
    g++ libatlas-base-dev gfortran libreadline-dev zlib1g-dev flex bison \
    libxml2-dev libxslt-dev libssl-dev -y
pip install virtualenv
apt-get build-dep python-lxml -y

# install sql servers
apt-get install postgresql-client libpq-dev -y

export PYTHONIOENCODING='utf-8'
export VIRTUALENV_DISTRIBUTE=0

function create_fake_pandas() {
    local site_pkg_dir="$1"
    rm -rf $site_pkg_dir/pandas
    mkdir $site_pkg_dir/pandas
    touch $site_pkg_dir/pandas/__init__.py
    echo "version = '0.10.0-phony'" > $site_pkg_dir/pandas/version.py
}


function get_site_pkgs_dir() {
    python$1 -c 'import distutils; print(distutils.sysconfig.get_python_lib())'
}


function create_wheel() {
    local pip_args="$1"
    local wheelhouse="$2"
    local n="$3"
    local pyver="$4"

    local site_pkgs_dir="$(get_site_pkgs_dir $pyver)"


    if [[ "$n" == *statsmodels* ]]; then
        create_fake_pandas $site_pkgs_dir && \
        pip wheel $pip_args --wheel-dir=$wheelhouse $n && \
        pip install $pip_args --no-index $n && \
        rm -Rf $site_pkgs_dir
    else
        pip wheel $pip_args --wheel-dir=$wheelhouse $n
        pip install $pip_args --no-index $n
    fi
}


function generate_wheels() {
    # get the requirements file
    local reqfile="$1"

    # get the python version
    local TAG=$(echo $reqfile |  grep -Po "(\d\.?[\d\-](_\w+)?)")

    # base dir for wheel dirs
    local WHEELSTREET=/wheelhouse
    local WHEELHOUSE="$WHEELSTREET/$TAG"

    local PY_VER="${TAG:0:3}"
    local PY_MAJOR="${PY_VER:0:1}"
    local PIP_ARGS="--use-wheel --find-links=$WHEELHOUSE --download-cache /tmp"

    # install the python version if not installed
    apt-get install python$PY_VER python$PY_VER-dev -y

    # create a new virtualenv
    rm -Rf /tmp/venv
    virtualenv -p python$PY_VER /tmp/venv
    source /tmp/venv/bin/activate

    # install pip setuptools
    pip install -I --download-cache /tmp 'git+https://github.com/pypa/pip@42102e9d#egg=pip'
    pip install -I -U --download-cache /tmp setuptools
    pip install -I --download-cache /tmp wheel

    # make the dir if it doesn't exist
    mkdir -p $WHEELHOUSE

    # put the requirements file in the wheelhouse
    cp $reqfile $WHEELHOUSE

    # install and build the wheels
    cat $reqfile | while read N; do
        create_wheel "$PIP_ARGS" "$WHEELHOUSE" "$N" "$PY_VER"
    done
}


# generate a single wheel version
# generate_wheels "/reqf/requirements-2.7.txt"
# 
# if vagrant is already up
# run as vagrant provision

for reqfile in $(ls -1 /reqf/requirements-*.*); do
    generate_wheels "$reqfile"
done
