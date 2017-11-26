#!/bin/bash

source activate pandas

echo "install 36 BUILD_TEST"

conda install -n pandas -c conda-forge pyarrow dask pyqt qtpy
