#!/bin/bash -e

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
