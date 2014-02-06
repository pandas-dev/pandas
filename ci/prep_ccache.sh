#!/bin/bash

if [ "$IRON_TOKEN" ]; then
    sudo apt-get $APT_ARGS install ccache p7zip-full
    # iron_cache, pending py3 fixes upstream
    pip install -I --allow-external --allow-insecure git+https://github.com/iron-io/iron_cache_python.git@8a451c7d7e4d16e0c3bedffd0f280d5d9bd4fe59#egg=iron_cache


    curdir=$(pwd)
    python ci/ironcache/get.py
    ccache -C

    if [ -f ~/ccache.7z ]; then
        echo "Cache retrieved"
        cd $HOME
        7za e $HOME/ccache.7z
        # ls -l $HOME
        cd /
        tar xvf $HOME/ccache
        rm -rf $HOME/ccache.7z
        rm -rf $HOME/ccache

    fi

    # did the last commit change cython files?
    cd $curdir

    echo "diff from HEAD~5"
    git diff HEAD~5 --numstat

    retval=$(git diff HEAD~5 --numstat | grep -P "pyx|pxd"|wc -l)
    echo "number of cython files changed: $retval"

    if [ $retval -eq 0 ]
    then
        # nope, reuse cython files
        echo "Will reuse cached cython file"
        touch "$TRAVIS_BUILD_DIR"/pandas/*.c
        touch "$TRAVIS_BUILD_DIR"/pandas/src/*.c
        touch "$TRAVIS_BUILD_DIR"/pandas/*.cpp
    else
        echo "Rebuilding cythonized files"
    fi
fi

exit 0
