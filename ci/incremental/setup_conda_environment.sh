#!/bin/bash

set -v -e

CONDA_INSTALL="conda install -q -y"
PIP_INSTALL="pip install -q"

# Deactivate any environment
source deactivate
# Display root environment (for debugging)
conda list
# Clean up any left-over from a previous build
# (note workaround for https://github.com/conda/conda/issues/2679:
#  `conda env remove` issue)
conda remove --all -q -y -n $CONDA_ENV

echo
echo "[create env]"
time conda env create -q -n "${CONDA_ENV}" --file="${ENV_FILE}" || exit 1

# Activate first
set +v
source activate $CONDA_ENV
set -v

# remove any installed pandas package
# w/o removing anything else
echo
echo "[removing installed pandas]"
conda remove pandas -y --force || true
pip uninstall -y pandas || true

echo
echo "[no installed pandas]"
conda list pandas

if [ -n "$LOCALE_OVERRIDE" ]; then
    sudo locale-gen "$LOCALE_OVERRIDE"
fi

# # Install the compiler toolchain
# if [[ $(uname) == Linux ]]; then
#     if [[ "$CONDA_SUBDIR" == "linux-32" || "$BITS32" == "yes" ]] ; then
#         $CONDA_INSTALL gcc_linux-32 gxx_linux-32
#     else
#         $CONDA_INSTALL gcc_linux-64 gxx_linux-64
#     fi
# elif  [[ $(uname) == Darwin ]]; then
#     $CONDA_INSTALL clang_osx-64 clangxx_osx-64
#     # Install llvm-openmp and intel-openmp on OSX too
#     $CONDA_INSTALL llvm-openmp intel-openmp
# fi
