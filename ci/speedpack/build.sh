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
    g++ libatlas-base-dev gfortran -y
pip install virtualenv
apt-get build-dep python-lxml -y

export PYTHONIOENCODING='utf-8'
export VIRTUALENV_DISTRIBUTE=0

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
        pip wheel $PIP_ARGS --wheel-dir=$WHEELHOUSE $N
        pip install $PIP_ARGS --no-index $N
    done
}


for reqfile in $(ls -1 /reqf/requirements-*.*); do
    generate_wheels "$reqfile"
done
