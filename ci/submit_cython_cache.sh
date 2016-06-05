#!/bin/bash

CACHE_File="$HOME/.cache/cython_files.tar"
home_dir=$(pwd)
ls "$TRAVIS_BUILD_DIR/pandas/"
tar cf "$CACHE_File" \
"$TRAVIS_BUILD_DIR"/pandas/{index,algos,lib,tslib,parser,hashtable}.c \
"$TRAVIS_BUILD_DIR"/pandas/src/{sparse,testing,period,period_helper}.c \
"$TRAVIS_BUILD_DIR"/pandas/msgpack/{_packer,_unpacker}.cpp

exit 0
