#!/bin/bash

home_dir=$(pwd)
ccache -s

MISSES=$(ccache -s | grep "cache miss" | grep -Po "\d+")
echo "MISSES: $MISSES"

if [ x"$MISSES" == x"0" ]; then
    echo "No cache misses detected, skipping upload"
    exit 0
fi

if [ "$IRON_TOKEN" ]; then

    # install the compiler cache
    sudo apt-get $APT_ARGS install ccache p7zip-full
    # iron_cache, pending py3 fixes upstream
    pip install -I --allow-external --allow-insecure git+https://github.com/iron-io/iron_cache_python.git@8a451c7d7e4d16e0c3bedffd0f280d5d9bd4fe59#egg=iron_cache

    rm -rf $HOME/ccache.7z

    tar cf - $HOME/.ccache \
    "$TRAVIS_BUILD_DIR"/pandas/{index,algos,lib,tslib,parser,hashtable}.c \
    "$TRAVIS_BUILD_DIR"/pandas/src/{sparse,testing}.c \
    "$TRAVIS_BUILD_DIR"/pandas/msgpack.cpp  \
    |  7za a -si $HOME/ccache.7z

    split -b 500000 -d $HOME/ccache.7z $HOME/ccache.

    python ci/ironcache/put.py
fi;

exit 0
