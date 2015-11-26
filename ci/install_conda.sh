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

function edit_init()
{
    if [ -n "$LOCALE_OVERRIDE" ]; then
        echo "Adding locale to the first line of pandas/__init__.py"
        rm -f pandas/__init__.pyc
        sedc="3iimport locale\nlocale.setlocale(locale.LC_ALL, '$LOCALE_OVERRIDE')\n"
        sed -i "$sedc" pandas/__init__.py
        echo "head -4 pandas/__init__.py"
        head -4 pandas/__init__.py
        echo
    fi
}

edit_init

python_major_version="${TRAVIS_PYTHON_VERSION:0:1}"
[ "$python_major_version" == "2" ] && python_major_version=""

home_dir=$(pwd)
echo "home_dir: [$home_dir]"

if [ -n "$LOCALE_OVERRIDE" ]; then
    # make sure the locale is available
    # probably useless, since you would need to relogin
    time sudo locale-gen "$LOCALE_OVERRIDE"
fi

# Need to enable for locale testing. The location of the locale file(s) is
# distro specific. For example, on Arch Linux all of the locales are in a
# commented file--/etc/locale.gen--that must be commented in to be used
# whereas Ubuntu looks in /var/lib/locales/supported.d/* and generates locales
# based on what's in the files in that folder
time echo 'it_CH.UTF-8 UTF-8' | sudo tee -a /var/lib/locales/supported.d/it
time sudo locale-gen


# install gui for clipboard testing
if [ -n "$CLIPBOARD_GUI" ]; then
    echo "Using CLIPBOARD_GUI: $CLIPBOARD_GUI"
    [ -n "$python_major_version" ] && py="py"
    python_cb_gui_pkg=python${python_major_version}-${py}${CLIPBOARD_GUI}
    time sudo apt-get $APT_ARGS install $python_cb_gui_pkg
fi


# install a clipboard if $CLIPBOARD is not empty
if [ -n "$CLIPBOARD" ]; then
    echo "Using clipboard: $CLIPBOARD"
    time sudo apt-get $APT_ARGS install $CLIPBOARD
fi

python_major_version="${TRAVIS_PYTHON_VERSION:0:1}"
[ "$python_major_version" == "2" ] && python_major_version=""

wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh || exit 1
bash miniconda.sh -b -p $HOME/miniconda || exit 1

conda config --set always_yes yes --set changeps1 no || exit 1
conda update -q conda || exit 1
conda config --add channels http://conda.anaconda.org/pandas || exit 1
conda config --set ssl_verify false || exit 1

# Useful for debugging any issues with conda
conda info -a || exit 1

# build deps
REQ="ci/requirements-${TRAVIS_PYTHON_VERSION}${JOB_TAG}.build"
time conda create -n pandas python=$TRAVIS_PYTHON_VERSION nose || exit 1
time conda install -n pandas --file=${REQ} || exit 1

source activate pandas

# set the compiler cache to work
if [ "$IRON_TOKEN" ]; then
    export PATH=/usr/lib/ccache:/usr/lib64/ccache:$PATH
    gcc=$(which gcc)
    echo "gcc: $gcc"
    ccache=$(which ccache)
    echo "ccache: $ccache"
    export CC='ccache gcc'
fi

if [ "$BUILD_TEST" ]; then

    # build testing
    pip uninstall --yes cython
    pip install cython==0.15.1
    ( python setup.py build_ext --inplace && python setup.py develop ) || true

else

    # build but don't install
    time python setup.py build_ext --inplace || exit 1

    # we may have run installations
    REQ="ci/requirements-${TRAVIS_PYTHON_VERSION}${JOB_TAG}.run"
    time conda install -n pandas --file=${REQ} || exit 1

    # we may have additional pip installs
    REQ="ci/requirements-${TRAVIS_PYTHON_VERSION}${JOB_TAG}.pip"
    if [ -e ${REQ} ]; then
        pip install -r $REQ
    fi

    # remove any installed pandas package
    conda remove pandas

    # install our pandas
    python setup.py develop  || exit 1

fi

true
