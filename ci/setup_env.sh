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

#!/bin/bash

set -v -e

# Install Miniconda
unamestr=`uname`
if [[ "$unamestr" == 'Linux' ]]; then
    if [[ "$BITS32" == "yes" ]]; then
        wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86.sh -O miniconda.sh
    else
        wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    fi
elif [[ "$unamestr" == 'Darwin' ]]; then
    wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh
else
  echo Error
fi
chmod +x miniconda.sh
./miniconda.sh -b


echo
echo "[show conda]"
which conda

echo
echo "[update conda]"
conda config --set ssl_verify false || exit 1
conda config --set quiet true --set always_yes true --set changeps1 false || exit 1
conda update -q conda

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

# Deactivate any environment
source deactivate
# Display root environment (for debugging)
conda list
# Clean up any left-over from a previous build
# (note workaround for https://github.com/conda/conda/issues/2679:
#  `conda env remove` issue)
conda remove --all -q -y -n pandas-dev

echo
echo "[create env]"

# create our environment
time conda env create -q --file="${ENV_FILE}" || exit 1

set +v
source activate pandas-dev
set -v

# remove any installed pandas package
# w/o removing anything else
echo
echo "[removing installed pandas]"
conda remove pandas -y --force || true
pip uninstall -y pandas || true

echo
echo "[no installed pandas]"
conda list pandas


if [ -n "$LOCALE_OVERRIDE" ]; then
    sudo locale-gen "$LOCALE_OVERRIDE"
fi

#!/bin/bash

# Make sure any error below is reported as such
set -v -e

echo "[building extensions]"
python setup.py build_ext -q --inplace
python -m pip install -e .

echo
echo "[show environment]"
conda list

echo
echo "[done]"
exit 0

#!/bin/bash

if [ "${TRAVIS_OS_NAME}" != "linux" ]; then
   echo "not using dbs on non-linux"
   exit 0
fi

echo "installing dbs"
mysql -e 'create database pandas_nosetest;'
psql -c 'create database pandas_nosetest;' -U postgres

echo "done"
exit 0
