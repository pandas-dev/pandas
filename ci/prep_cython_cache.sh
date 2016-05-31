#!/bin/bash

CACHE_File="$HOME/.cache/cython_files.tar"

clear_cache=0

if [ -f "$CACHE_File" ]; then

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

if [ $clear_cache -eq 1 ] && [ $retval -eq 0 ]
then
    # nope, reuse cython files
    echo "Will reuse cached cython file"
    touch "$TRAVIS_BUILD_DIR"/pandas/*.c
    touch "$TRAVIS_BUILD_DIR"/pandas/src/*.c
    touch "$TRAVIS_BUILD_DIR"/pandas/*.cpp
else
    echo "Rebuilding cythonized files"
fi


exit 0
