#!/bin/bash

source activate $CONDA_ENV

# Make sure any error below is reported as such
set -v -e

echo "[building extensions]"
python setup.py build_ext -q --inplace
python -m pip install -e

echo
echo "[show environment]"
conda list

echo
echo "[done]"
exit 0
