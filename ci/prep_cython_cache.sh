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
    if [ "$TRAVIS_PULL_REQUEST" == "false" ]
    then
        echo "Not a PR: checking for cython files changes from last 2 commits"
        git diff HEAD~2 --numstat | grep -E "pyx|pxd"
        retval=$(git diff HEAD~2 --numstat | grep -E "pyx|pxd"| wc -l)
    else
        echo "PR: checking for any cython file changes from last 5 commits"
        git diff PR_HEAD~5 --numstat | grep -E "pyx|pxd"
        retval=$(git diff PR_HEAD~5 --numstat | grep -E "pyx|pxd"| wc -l)
    fi
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
