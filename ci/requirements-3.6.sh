#!/bin/bash

source activate pandas

echo "install 36"

conda install -n pandas -c conda-forge feather-format
