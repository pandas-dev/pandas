#!/bin/bash

source activate pandas

echo "install 36 CONDA_BUILD_TEST"

conda install conda-build

conda build ./conda.recipe/ --numpy 1.11 --python 3.6 -q
