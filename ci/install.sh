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

# known working
# pip==1.5.1
# setuptools==2.2
# wheel==0.22
# nose==1.3.0 (1.3.1 broken for PY3)

pip install -I -U pip
pip install -I -U setuptools
pip install wheel==0.22

# install nose
pip uninstall nose -y

if [ -n "$EXPERIMENTAL" ]; then

    # install from master
    rm -Rf /tmp/nose
    cd /tmp
    git clone --branch master https://github.com/nose-devs/nose.git nose
    cd nose
    python setup.py install
    cd $home_dir

else

    # known good version
    pip install nose==1.3.0

fi


# comment this line to disable the fetching of wheel files
base_url=http://pandas.pydata.org/pandas-build/dev/wheels

wheel_box=${TRAVIS_PYTHON_VERSION}${JOB_TAG}
PIP_ARGS+=" -I --use-wheel --find-links=$base_url/$wheel_box/ --allow-external --allow-insecure"

if [ -n "$LOCALE_OVERRIDE" ]; then
    # make sure the locale is available
    # probably useless, since you would need to relogin
    time sudo locale-gen "$LOCALE_OVERRIDE"
fi

# we need these for numpy
time sudo apt-get $APT_ARGS install libatlas-base-dev gfortran

if [ -n "$NUMPY_BUILD" ]; then
    # building numpy

    cd $home_dir
    echo "cloning numpy"

    rm -Rf /tmp/numpy
    cd /tmp

    # remove the system installed numpy
    pip uninstall numpy -y

    # install cython
    pip install --find-links http://wheels.astropy.org/ --find-links http://wheels2.astropy.org/ --use-wheel Cython

    # clone & install
    git clone --branch $NUMPY_BUILD https://github.com/numpy/numpy.git numpy
    cd numpy
    time pip install .
    pip uninstall cython -y

    cd $home_dir
    numpy_version=$(python -c 'import numpy; print(numpy.__version__)')
    echo "[$home_dir] numpy current: $numpy_version"
fi

# Force virtualenv to accept system_site_packages
rm -f $VIRTUAL_ENV/lib/python$TRAVIS_PYTHON_VERSION/no-global-site-packages.txt

time pip install $PIP_ARGS -r ci/requirements-${wheel_box}.txt

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


# Optional Deps
if [ -n "$FULL_DEPS" ]; then
    echo "Installing FULL_DEPS"

    # need libhdf5 for PyTables
    time sudo apt-get $APT_ARGS install libhdf5-serial-dev
fi


# set the compiler cache to work
if [ "$IRON_TOKEN" ]; then
    export PATH=/usr/lib/ccache:/usr/lib64/ccache:$PATH
    gcc=$(which gcc)
    echo "gcc: $gcc"
    ccache=$(which ccache)
    echo "ccache: $ccache"
    export CC='ccache gcc'
fi

# build pandas
time python setup.py sdist
pip uninstall cython -y

# install pandas
time pip install $(find dist | grep gz | head -n 1)

# restore cython (if not numpy building)
if [ -z "$NUMPY_BUILD" ]; then
    time pip install $PIP_ARGS  $(cat ci/requirements-${wheel_box}.txt | grep -i cython)
fi

true
