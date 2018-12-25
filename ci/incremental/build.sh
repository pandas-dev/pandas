#!/bin/bash

# Make sure any error below is reported as such
set -v -e

echo "[building extensions]"
python setup.py build_ext -q --inplace
python -m pip install -e .

echo -e "\n[show environment]"
conda list

echo -e "\n[done]"
exit 0
