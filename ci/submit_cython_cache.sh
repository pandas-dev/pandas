#!/bin/bash

CACHE_File="$HOME/.cache/cython_files.tar"
home_dir=$(pwd)
ls "$TRAVIS_BUILD_DIR"/pandas/"
tar cf $CACHE_DIR/cython_files.tar \
"$TRAVIS_BUILD_DIR"/pandas/{index,algos,lib,tslib,parser,hashtable}.c \
"$TRAVIS_BUILD_DIR"/pandas/src/{sparse,testing,period,}.c \
"$TRAVIS_BUILD_DIR"/pandas/msgpack/{_packer,_unpacker,period_helper}.cpp"

exit 0
