#!/bin/bash

# edit the locale file if needed
function edit_init()
{
    if [ -n "$LOCALE_OVERRIDE" ]; then
        echo "[Adding locale to the first line of pandas/__init__.py]"
        rm -f pandas/__init__.pyc
        sedc="3iimport locale\nlocale.setlocale(locale.LC_ALL, '$LOCALE_OVERRIDE')\n"
        sed -i "$sedc" pandas/__init__.py
        echo "[head -4 pandas/__init__.py]"
        head -4 pandas/__init__.py
        echo
    fi
}

echo
echo "[install_travis]"
edit_init

home_dir=$(pwd)
echo
echo "[home_dir]: $home_dir"

# install miniconda
MINICONDA_DIR="$HOME/miniconda3"

echo
echo "[Using clean Miniconda install]"

if [ -d "$MINICONDA_DIR" ]; then
    rm -rf "$MINICONDA_DIR"
fi

# install miniconda
if [ "${TRAVIS_OS_NAME}" == "osx" ]; then
    time wget http://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh || exit 1
else
    time wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh || exit 1
fi
time bash miniconda.sh -b -p "$MINICONDA_DIR" || exit 1

echo
echo "[show conda]"
which conda

echo
echo "[update conda]"
conda config --set ssl_verify false || exit 1
conda config --set always_yes true --set changeps1 false || exit 1
conda update -q conda

echo
echo "[add channels]"
conda config --remove channels defaults || exit 1
conda config --add channels defaults || exit 1

if [ "$CONDA_FORGE" ]; then
    # add conda-forge channel as priority
    conda config --add channels conda-forge || exit 1
fi

# Useful for debugging any issues with conda
conda info -a || exit 1

# set the compiler cache to work
echo
if [ -z "$NOCACHE" ] && [ "${TRAVIS_OS_NAME}" == "linux" ]; then
    echo "[Using ccache]"
    export PATH=/usr/lib/ccache:/usr/lib64/ccache:$PATH
    gcc=$(which gcc)
    echo "[gcc]: $gcc"
    ccache=$(which ccache)
    echo "[ccache]: $ccache"
    export CC='ccache gcc'
elif [ -z "$NOCACHE" ] && [ "${TRAVIS_OS_NAME}" == "osx" ]; then
    echo "[Install ccache]"
    brew install ccache > /dev/null 2>&1
    echo "[Using ccache]"
    export PATH=/usr/local/opt/ccache/libexec:$PATH
    gcc=$(which gcc)
    echo "[gcc]: $gcc"
    ccache=$(which ccache)
    echo "[ccache]: $ccache"
else
    echo "[Not using ccache]"
fi

echo
echo "[create env]"

# create our environment
REQ="ci/requirements-${JOB}.build"
time conda create -n pandas --file=${REQ} || exit 1

source activate pandas

# may have addtl installation instructions for this build
echo
echo "[build addtl installs]"
REQ="ci/requirements-${JOB}.build.sh"
if [ -e ${REQ} ]; then
    time bash $REQ || exit 1
fi

time conda install -n pandas pytest
time pip install pytest-xdist

if [ "$LINT" ]; then
   conda install flake8
   pip install cpplint
fi

if [ "$COVERAGE" ]; then
    pip install coverage pytest-cov
fi

echo
if [ -z "$BUILD_TEST" ]; then

    # build but don't install
    echo "[build em]"
    time python setup.py build_ext --inplace || exit 1

fi

# we may have run installations
echo
echo "[conda installs]"
REQ="ci/requirements-${JOB}.run"
if [ -e ${REQ} ]; then
    time conda install -n pandas --file=${REQ} || exit 1
fi

# we may have additional pip installs
echo
echo "[pip installs]"
REQ="ci/requirements-${JOB}.pip"
if [ -e ${REQ} ]; then
   pip install -r $REQ
fi

# may have addtl installation instructions for this build
echo
echo "[addtl installs]"
REQ="ci/requirements-${JOB}.sh"
if [ -e ${REQ} ]; then
    time bash $REQ || exit 1
fi

# remove any installed pandas package
# w/o removing anything else
echo
echo "[removing installed pandas]"
conda remove pandas -y --force
pip uninstall -y pandas

if [ "$BUILD_TEST" ]; then

    # remove any installation
    pip uninstall -y pandas
    conda list pandas
    pip list --format columns |grep pandas

    # build & install testing
    echo ["building release"]
    bash scripts/build_dist_for_release.sh
    conda uninstall -y cython
    time pip install dist/*tar.gz || exit 1

else

    # install our pandas
    echo
    echo "[running setup.py develop]"
    python setup.py develop  || exit 1

fi

echo
echo "[show pandas]"
conda list pandas

echo
echo "[done]"
exit 0
