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
conda install pip conda  # create conda to create a historical artifact for pip & setuptools
conda update -n base conda

echo "conda info -a"
conda info -a

echo "source deactivate"
source deactivate

echo "conda list (root environment)"
conda list

# Clean up any left-over from a previous build
conda remove --all -q -y -n pandas-dev

echo
echo "conda env create -q --file=${ENV_FILE}"
time conda env create -q --file="${ENV_FILE}"


if [[ "$BITS32" == "yes" ]]; then
    # activate 32-bit compiler
    export CONDA_BUILD=1
fi

echo "activate pandas-dev"
source activate pandas-dev

# Explicitly set an environment variable indicating that this is pandas' CI environment.
#
# This allows us to enable things like -Werror that shouldn't be activated in
# downstream CI jobs that may also build pandas from source.
export PANDAS_CI=1

echo
echo "remove any installed pandas package"
echo "w/o removing anything else"
conda remove pandas -y --force || true
pip uninstall -y pandas || true

echo
echo "remove postgres if has been installed with conda"
echo "we use the one from the CI"
conda remove postgresql -y --force || true

echo
echo "remove qt"
echo "causes problems with the clipboard, we use xsel for that"
conda remove qt -y --force || true

echo
echo "conda list pandas"
conda list pandas

# Make sure any error below is reported as such

echo "[Build extensions]"
python setup.py build_ext -q -j2

echo "[Updating pip]"
python -m pip install --no-deps -U pip wheel setuptools

echo "[Install pandas]"
python -m pip install --no-build-isolation -e .

echo
echo "conda list"
conda list

# Install DB for Linux

if [[ -n ${SQL:0} ]]; then
  echo "installing dbs"
  mysql -e 'create database pandas_nosetest;'
  psql -c 'create database pandas_nosetest;' -U postgres
else
   echo "not using dbs on non-linux Travis builds or Azure Pipelines"
fi
echo "done"
