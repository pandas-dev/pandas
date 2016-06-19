#!/bin/bash

ls "$HOME/.cache/"
CACHE_File="$HOME/.cache/cython_files.tar"

clear_cache=0
home_dir=$(pwd)

if [ -f "$CACHE_File" ] && [ "$USE_CACHE" ]; then

    echo "Cache available"
    clear_cache=1
    # did the last commit change cython files?
    # go back 2 commits
    retval=$(git diff HEAD~2 --numstat | grep -E "pyx|pxd"| wc -l)
    echo "number of cython files changed: $retval"
fi

if [ $clear_cache -eq 1 ] && [ $retval -eq 0 ] && [ "$USE_CACHE" ]
then
    # nope, reuse cython files
    echo "Will reuse cached cython file"
    cd /
    tar xvmf $CACHE_File
    cd $home_dir
else
    echo "Rebuilding cythonized files"
    echo "Use cache = $USE_CACHE"
    echo "Clear cache = $clear_cache"
fi


exit 0
