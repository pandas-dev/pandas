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


MINICONDA_DIR=/usr/local/miniconda
if [ -e $MINICONDA_DIR ] && [ "$BITS32" != yes ]; then
    echo "Found Miniconda installation at $MINICONDA_DIR"
else
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
fi
export PATH=$MINICONDA_DIR/bin:$PATH

echo
echo "which conda"
which conda

echo
echo "update conda"
conda config --set ssl_verify false
conda config --set quiet true --set always_yes true --set changeps1 false
conda install -y -c conda-forge -n base 'mamba>=0.21.2' pip setuptools

echo "conda info -a"
conda info -a

echo "conda list (root environment)"
conda list

echo
# Clean up any left-over from a previous build
mamba env remove -n pandas-dev
echo "mamba env update --file=${ENV_FILE}"
# See https://github.com/mamba-org/mamba/issues/633
mamba create -q -n pandas-dev
time mamba env update -n pandas-dev --file="${ENV_FILE}"

echo "conda list -n pandas-dev"
conda list -n pandas-dev

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

if pip list | grep -q ^pandas; then
    echo
    echo "remove any installed pandas package w/o removing anything else"
    pip uninstall -y pandas || true
fi

if [ "$(conda list -f qt --json)" != [] ]; then
    echo
    echo "remove qt"
    echo "causes problems with the clipboard, we use xsel for that"
    conda remove qt -y --force || true
fi

echo "Build extensions"
python setup.py build_ext -q -j3

echo "Install pandas"
python -m pip install --no-build-isolation --no-use-pep517 -e .

echo "done"
