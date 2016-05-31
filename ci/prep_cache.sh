#!/bin/bash

CACHE_DIR="$HOME/.cache"

clear_cache=0

if [ -f "$CACHE_DIR/cache.tar.gz" ]; then

    home_dir=$(pwd)

    echo "Cache retrieved"
    clear_cache=1
    cd $HOME
    # ls -l $HOME
    cd /
    tar xzvf $CACHE_DIR/cache.tar.gz

    # did the last commit change cython files?
    cd $home_dir

    retval=$(git diff HEAD~3 --numstat | grep -P "pyx|pxd"|wc -l)
    echo "number of cython files changed: $retval"

fi

if [ $clear_cache -eq 1 ] && [ $retval -eq 0 ]
then
    # nope, reuse cython files
    echo "Will reuse cached cython file"
    touch "$TRAVIS_BUILD_DIR"/pandas/*.c
    touch "$TRAVIS_BUILD_DIR"/pandas/src/*.c
    touch "$TRAVIS_BUILD_DIR"/pandas/*.cpp
else
    echo "Rebuilding cythonized files"
    rm -rf $CACHE_DIR/cache.tar.gz
fi


exit 0
