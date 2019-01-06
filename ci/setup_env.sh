#!/bin/bash -e


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
        sudo locale-gen "$LOCALE_OVERRIDE"
    fi
}

edit_init

home_dir=$(pwd)
echo
echo "home_dir: $home_dir"

# install miniconda
MINICONDA_DIR="$HOME/miniconda3"

echo
echo "Using clean Miniconda install"

if [ -d "$MINICONDA_DIR" ]; then
    rm -rf "$MINICONDA_DIR"
fi

echo "Install Miniconda"
unamestr=`uname`
if [[ "$unamestr" == 'Linux' ]]; then
    if [[ "$BITS32" == "yes" ]]; then
        conda_dist="Linux-x86.sh"
    else
        conda_dist="Linux-x86_64.sh"
    fi
elif [[ "$unamestr" == 'Darwin' ]]; then
    conda_dist="MacOSX-x86_64.sh"
else
  echo "OS $unamestr not supported"
  exit 1
fi

wget -q "https://repo.continuum.io/miniconda/Miniconda3-latest-$conda_dist" -O miniconda.sh
chmod +x miniconda.sh
./miniconda.sh -b

echo
echo "Show conda"
which conda

echo
echo "Update conda"
conda config --set ssl_verify false
conda config --set quiet true --set always_yes true --set changeps1 false
conda update -q conda

# Useful for debugging any issues with conda
conda info -a

echo
echo "set the compiler cache to work"
if [ -z "$NOCACHE" ] && [ "${TRAVIS_OS_NAME}" == "linux" ]; then
    echo "Using ccache"
    export PATH=/usr/lib/ccache:/usr/lib64/ccache:$PATH
    gcc=$(which gcc)
    echo "gcc: $gcc"
    ccache=$(which ccache)
    echo "ccache: $ccache"
    export CC='ccache gcc'
elif [ -z "$NOCACHE" ] && [ "${TRAVIS_OS_NAME}" == "osx" ]; then
    echo "Install ccache"
    brew install ccache > /dev/null 2>&1
    echo "Using ccache"
    export PATH=/usr/local/opt/ccache/libexec:$PATH
    gcc=$(which gcc)
    echo "gcc: $gcc"
    ccache=$(which ccache)
    echo "ccache: $ccache"
else
    echo "Not using ccache"
fi

echo "Deactivate any environment"
source deactivate

echo "Display root environment (for debugging)"
conda list

# Clean up any left-over from a previous build
# (note workaround for https://github.com/conda/conda/issues/2679:
#  `conda env remove` issue)
conda remove --all -q -y -n pandas-dev

# create our environment
echo
echo "Create environment"
time conda env create -q --file="${ENV_FILE}"

set +v
source activate pandas-dev
set -v

# remove any installed pandas package
# w/o removing anything else
echo
echo "Removing installed pandas"
conda remove pandas -y --force || true
pip uninstall -y pandas || true

echo
echo "No installed pandas"
conda list pandas

# Make sure any error below is reported as such

echo "Build extensions and install pandas"
python setup.py build_ext -q --inplace
python -m pip install -e .

echo
echo "Show environment"
conda list

# Install DB for Linux
if [ "${TRAVIS_OS_NAME}" != "linux" ]; then
   echo "not using dbs on non-linux"
else
   echo "installing dbs"
   mysql -e 'create database pandas_nosetest;'
   psql -c 'create database pandas_nosetest;' -U postgres
fi

echo
echo "Setup virtual framebuffer"
export DISPLAY=":99.
if [ "${TRAVIS_OS_NAME}" == "linux" ]; then
   sh -e /etc/init.d/xvfb start
   sleep 3
fi

echo "done"
exit 0
