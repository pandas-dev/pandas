#!/usr/bin/env bash

home_dir=$(pwd)
echo "[home_dir: $home_dir]"

echo "[ls -ltr]"
ls -ltr

echo "[Using clean Miniconda install]"
rm -rf "$MINICONDA_DIR"

# install miniconda
wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -q -O miniconda.sh || exit 1
bash miniconda.sh -b -p "$MINICONDA_DIR" || exit 1

export PATH="$MINICONDA_DIR/bin:$PATH"

echo "[update conda]"
conda config --set ssl_verify false || exit 1
conda config --set always_yes true --set changeps1 false || exit 1
conda update -q conda

# add the pandas channel to take priority
# to add extra packages
echo "[add channels]"
conda config --add channels pandas || exit 1
conda config --remove channels defaults || exit 1
conda config --add channels defaults || exit 1

# Useful for debugging any issues with conda
conda info -a || exit 1

# support env variables passed
export ENVS_FILE=".envs"

# make sure that the .envs file exists. it is ok if it is empty
touch $ENVS_FILE

# assume all command line arguments are environmental variables
for var in "$@"
do
    echo "export $var" >> $ENVS_FILE
done

echo "[environmental variable file]"
cat $ENVS_FILE
source $ENVS_FILE

# edit the locale override if needed
if [ -n "$LOCALE_OVERRIDE" ]; then
    echo "[Adding locale to the first line of pandas/__init__.py]"
    rm -f pandas/__init__.pyc
    sedc="3iimport locale\nlocale.setlocale(locale.LC_ALL, '$LOCALE_OVERRIDE')\n"
    sed -i "$sedc" pandas/__init__.py
    echo "[head -4 pandas/__init__.py]"
    head -4 pandas/__init__.py
    echo
fi

# create envbuild deps
echo "[create env]"
time conda env create -q -n pandas --file="${ENV_FILE}" || exit 1

source activate pandas

# remove any installed pandas package
# w/o removing anything else
echo
echo "[removing installed pandas]"
conda remove pandas -y --force
pip uninstall -y pandas

# build but don't install
echo "[build em]"
time python setup.py build_ext --inplace || exit 1

echo
echo "[show environment]"

conda list
