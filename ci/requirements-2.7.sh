#!/bin/bash

source activate pandas

echo "install 27"

conda install -n pandas -c conda-forge feather-format pyarrow=0.4.1 jemalloc=4.5.0.post fastparquet
