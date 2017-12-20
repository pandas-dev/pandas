#!/bin/bash

source activate pandas

echo "install 36 PIP_BUILD_TEST"

conda install -n pandas -c conda-forge pyarrow dask pyqt qtpy
