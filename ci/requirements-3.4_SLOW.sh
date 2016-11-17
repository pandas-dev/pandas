#!/bin/bash

source activate pandas

echo "install 34_slow"

conda install -n pandas -c conda-forge/label/rc -c conda-forge matplotlib
