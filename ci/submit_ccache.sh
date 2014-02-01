#!/bin/bash

ccache -s

if [ "$IRON_TOKEN" ]; then
    rm -rf ~/ccache.7z

    tar cf - $HOME/.ccache \
        |  7za a -si ~/ccache.7z

    # Caching the cython files could be be made to work.
    # There's a timestamp issue to be handled (http://stackoverflow.com/questions/1964470)
    # but you could get the build to reuse them.
    # However, there's a race condition between travis cythonizing
    # and a new commit being made on a GH server somewhere.
    # and there's convenient causal link between the two we can
    # take advantage of.
    # Let's just wait another 60 seconds.
    # "$TRAVIS_BUILD_DIR"/pandas/{index,algos,lib,tslib,parser,hashtable}.c \
    # "$TRAVIS_BUILD_DIR"/pandas/src/{sparse,testing}.c \
    # "$TRAVIS_BUILD_DIR"/pandas/msgpack.cpp  \
    # |  7za a -si ~/ccache.7z

    split -b 500000 -d ~/ccache.7z ~/ccache.

    python ci/ironcache/put.py
fi;

exit 0
