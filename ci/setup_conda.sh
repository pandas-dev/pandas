#!/bin/bash -e

# edit the locale file if needed
if [[ "$(uname)" == "Linux" && -n "$LC_ALL" ]]; then
    echo "Adding locale to the first line of pandas/__init__.py"
    rm -f pandas/__init__.pyc
    SEDC="3iimport locale\nlocale.setlocale(locale.LC_ALL, '$LC_ALL')\n"
    sed -i "$SEDC" pandas/__init__.py

    echo "[head -4 pandas/__init__.py]"
    head -4 pandas/__init__.py
    echo
fi

echo "Install Miniconda"
DEFAULT_CONDA_URL="https://repo.continuum.io/miniconda/Miniconda3-latest"
if [[ "$(uname -m)" == 'aarch64' ]]; then
    CONDA_URL="https://github.com/conda-forge/miniforge/releases/download/4.10.1-4/Miniforge3-4.10.1-4-Linux-aarch64.sh"
elif [[ "$(uname)" == 'Linux' ]]; then
    if [[ "$BITS32" == "yes" ]]; then
        CONDA_URL="$DEFAULT_CONDA_URL-Linux-x86.sh"
    else
        CONDA_URL="$DEFAULT_CONDA_URL-Linux-x86_64.sh"
    fi
elif [[ "$(uname)" == 'Darwin' ]]; then
    CONDA_URL="$DEFAULT_CONDA_URL-MacOSX-x86_64.sh"
else
    echo "OS $(uname) not supported"
    exit 1
fi
echo "Downloading $CONDA_URL"
wget -q $CONDA_URL -O miniconda.sh
chmod +x miniconda.sh

MINICONDA_DIR="$HOME/miniconda3"
rm -rf $MINICONDA_DIR
./miniconda.sh -b -p $MINICONDA_DIR
export PATH=$MINICONDA_DIR/bin:$PATH

echo
echo "which conda"
which conda

echo
echo "update conda"
conda config --set ssl_verify false
conda config --set quiet true --set always_yes true --set changeps1 false
conda install pip conda # create conda to create a historical artifact for pip & setuptools
conda update -n base conda

echo "conda info -a"
conda info -a

echo "source deactivate"
source deactivate

echo "conda list (root environment)"
conda list

# Clean up any left-over from a previous build
conda remove --all -q -y -n pandas-dev
