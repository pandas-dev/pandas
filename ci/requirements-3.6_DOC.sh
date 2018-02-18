#!/bin/bash

source activate pandas

echo "[install DOC_BUILD deps]"

pip install pandas-gbq

conda install -n pandas -c conda-forge sphinx feather-format pyarrow nbsphinx pandoc fastparquet

conda install -n pandas -c r r rpy2 --yes
