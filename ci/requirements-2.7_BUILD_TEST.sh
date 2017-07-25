#!/bin/bash

source activate pandas

echo "install 27 BUILD_TEST"

conda install -n pandas -c conda-forge pyarrow=0.4.1 dask
