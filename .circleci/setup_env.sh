#!/bin/bash -e

echo "Install Mambaforge"
MAMBA_URL="https://github.com/conda-forge/miniforge/releases/download/4.14.0-0/Mambaforge-4.14.0-0-Linux-aarch64.sh"
echo "Downloading $MAMBA_URL"
wget -q $MAMBA_URL -O minimamba.sh
chmod +x minimamba.sh

MAMBA_DIR="$HOME/miniconda3"
rm -rf $MAMBA_DIR
./minimamba.sh -b -p $MAMBA_DIR

export PATH=$MAMBA_DIR/bin:$PATH

echo
echo "which conda"
which conda

echo
echo "update conda"
conda config --set ssl_verify false
conda config --set quiet true --set always_yes true --set changeps1 false
mamba install -y -c conda-forge -n base pip setuptools

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

echo "Build extensions"
# GH 47305: Parallel build can causes flaky ImportError from pandas/_libs/tslibs
python setup.py build_ext -q -j1

echo "Install pandas"
python -m pip install --no-build-isolation --no-use-pep517 -e .

echo "done"
