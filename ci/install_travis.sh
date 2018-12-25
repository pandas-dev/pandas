#!/bin/bash

# edit the locale file if needed
function edit_init()
{
    if [ -n "$LOCALE_OVERRIDE" ]; then
        echo "[Adding locale to the first line of pandas/__init__.py]"
        rm -f pandas/__init__.pyc
        sedc="3iimport locale\nlocale.setlocale(locale.LC_ALL, '$LOCALE_OVERRIDE')\n"
        sed -i "$sedc" pandas/__init__.py
        echo -e "[head -4 pandas/__init__.py]\n"
        head -4 pandas/__init__.py
    fi
}

echo -e "\n[install_travis]"
edit_init

home_dir=$(pwd)
echo -e "\n[home_dir]: $home_dir"

# install miniconda
MINICONDA_DIR="${HOME}/miniconda3"

echo -e "\n[Using clean Miniconda install]"

if [ -d "${MINICONDA_DIR}" ]; then
    rm -rf "${MINICONDA_DIR}"
fi

# install miniconda
if [ "${TRAVIS_OS_NAME}" == "osx" ]; then
    time wget http://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -q -O miniconda.sh || exit 1
else
    time wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -q -O miniconda.sh || exit 1
fi
time bash miniconda.sh -b -p "${MINICONDA_DIR}" || exit 1

echo -e "\n[show conda]"
which conda

echo -e "\n[update conda]"
conda config --set ssl_verify false || exit 1
conda config --set quiet true --set always_yes true --set changeps1 false || exit 1
conda update -q conda

# Useful for debugging any issues with conda
conda info -a || exit 1

# set the compiler cache to work
if [ -z "${NOCACHE}" ] && [ "${TRAVIS_OS_NAME}" == "linux" ]; then
    echo -e "\n[Using ccache]"
    export PATH=/usr/lib/ccache:/usr/lib64/ccache:$PATH
    echo "[gcc]: $(command -v gcc)"
    echo "[ccache]: $(command -v ccache)"
    export CC='ccache gcc'
elif [ -z "${NOCACHE}" ] && [ "${TRAVIS_OS_NAME}" == "osx" ]; then
    echo "[Install ccache]"
    brew install ccache > /dev/null 2>&1
    echo "[Using ccache]"
    export PATH=/usr/local/opt/ccache/libexec:$PATH
    echo "[gcc]: $(command -v gcc)"
    echo "[ccache]: $(command -v ccache)"
else
    echo "[Not using ccache]"
fi

echo -e "\n[create env]"

# create our environment
time conda env create -q --file="${ENV_FILE}" || exit 1

source activate pandas-dev

# remove any installed pandas package
# w/o removing anything else
echo -e "\n[removing installed pandas]"
conda remove pandas -y --force
pip uninstall -y pandas

echo -e "\n[no installed pandas]"
conda list pandas
pip list --format columns | grep pandas

# build and install
echo "[running setup.py develop]"
python setup.py develop  || exit 1

echo -e "\\n[show environment]"
conda list

echo -e "\n[done]"
exit 0
