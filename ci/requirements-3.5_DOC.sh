#!/bin/bash

source activate pandas

echo "[install DOC_BUILD deps]"

pip install pandas-gbq

conda install -n pandas -c conda-forge feather-format

conda install -n pandas -c r r rpy2 --yes
