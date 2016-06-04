#!/bin/bash

ls "$HOME/.cache/"
CACHE_File="$HOME/.cache/cython_files.tar"

clear_cache=0

if [ -f "$CACHE_File" ] && [ "$USE_CACHE" ]; then

    home_dir=$(pwd)

    echo "Cache retrieved"
    clear_cache=1
    cd $HOME
    # ls -l $HOME
    cd /
    tar xvf $CACHE_File

    # did the last commit change cython files?
    cd $home_dir

    retval=$(git diff HEAD~3 --numstat | grep -P "pyx|pxd"|wc -l)
    echo "number of cython files changed: $retval"

    rm -rf $CACHE_File

fi

if [ $clear_cache -eq 1 ] && [ $retval -eq 0 ] && [ "$USE_CACHE" ]
then
    # nope, reuse cython files
    echo "Will reuse cached cython file"
    touch "$TRAVIS_BUILD_DIR"/pandas/*.c
    touch "$TRAVIS_BUILD_DIR"/pandas/src/*.c
    touch "$TRAVIS_BUILD_DIR"/pandas/*.cpp
    touch "$TRAVIS_BUILD_DIR"/pandas/msgpack/*.cpp
else
    echo "Rebuilding cythonized files"
    echo "Use cache = $USE_CACHE"
    echo "Clear cache = $clear_cache"
fi


exit 0
