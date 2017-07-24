#!/bin/bash

source activate pandas

echo "[install DOC_BUILD deps]"

pip install pandas-gbq

conda install -n pandas -c conda-forge feather-format pyarrow=0.4.1 nbsphinx pandoc

conda install -n pandas -c r r rpy2 --yes
