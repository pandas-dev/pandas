#!/bin/bash

source activate pandas

echo "[install DOC_BUILD deps]"

pip install pandas-gbq

conda install -n pandas -c conda-forge feather-format pyarrow nbsphinx pandoc fastparquet=0.1.3

conda install -n pandas -c r r rpy2 --yes
