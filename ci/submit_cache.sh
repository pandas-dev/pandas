#!/bin/bash

CACHE_DIR="$HOME/.cache"
home_dir=$(pwd)

if [ -f "$CACHE_DIR/cache.tar.gz" ]; then

    echo "Not creating cache - file exists"

else

    ls "$TRAVIS_BUILD_DIR"/pandas/"
    tar cfz $CACHE_DIR/cache.tar.gz \
    "$TRAVIS_BUILD_DIR"/pandas/{index,algos,lib,tslib,parser,hashtable}.c \
    "$TRAVIS_BUILD_DIR"/pandas/src/{sparse,testing}.c \
    "$TRAVIS_BUILD_DIR"/pandas/msgpack.cpp

fi

exit 0
