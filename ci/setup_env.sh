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


echo "Install Mambaforge"
DEFAULT_MAMBA_URL="https://github.com/conda-forge/miniforge/releases/latest/download"
if [[ "$(uname -m)" == 'aarch64' ]]; then
    MAMBA_NAME="Mambaforge-Linux-aarch64.sh"
    MAMBA_URL="$DEFAULT_MAMBA_URL/$MAMBA_NAME"
elif [[ "$(uname)" == 'Linux' ]]; then
    if [[ "$BITS32" == "yes" ]]; then
        MAMBA_NAME="Miniconda3-latest-Linux-x86.sh"
        MAMBA_URL="https://repo.continuum.io/miniconda/$MAMBA_NAME"
    else
        MAMBA_NAME="Mambaforge-Linux-x86_64.sh"
        MAMBA_URL="$DEFAULT_MAMBA_URL/$MAMBA_NAME"
    fi
elif [[ "$(uname)" == 'Darwin' ]]; then
    MAMBA_NAME="Mambaforge-MacOSX-x86_64.sh"
    MAMBA_URL="$DEFAULT_MAMBA_URL/$MAMBA_NAME"
else
  echo "OS $(uname) not supported"
  exit 1
fi
echo "Downloading $MAMBA_URL"
wget -q $MAMBA_URL -O $MAMBA_NAME
chmod +x $MAMBA_NAME

MAMBA_DIR="$HOME/mambaforge"
echo "Mamba directory $MAMBA_DIR"
rm -rf $MAMBA_DIR
./$MAMBA_NAME -b -p $MAMBA_DIR
export PATH=$MAMBA_DIR/bin:$PATH
export PATH=$MAMBA_DIR/condabin:$PATH

echo
echo "which mamba"
which mamba

echo
echo "update mamba"
conda config --set ssl_verify false
conda config --set quiet true --set always_yes true --set changeps1 false
mamba install -y pip conda  # create mamba to create a historical artifact for pip & setuptools
mamba update -y -n base conda

echo "mamba info -a"
mamba info -a

echo "source deactivate"
source deactivate

echo "mamba list (root environment)"
mamba list

# Clean up any left-over from a previous build
conda remove --all -q -y -n pandas-dev

echo
echo "mamba env create -q --file=${ENV_FILE}"
time mamba env create -q --file="${ENV_FILE}"


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
