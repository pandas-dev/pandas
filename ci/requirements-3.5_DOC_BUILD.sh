#!/bin/bash

source activate pandas

echo "install DOC_BUILD"

conda install -n pandas -c conda-forge feather-format
