#!/bin/bash

ls "$HOME/.cache/"

PYX_CACHE_DIR="$HOME/.cache/pyxfiles"
pyx_file_list=`find ${TRAVIS_BUILD_DIR} -name "*.pyx"`
CACHE_File="$HOME/.cache/cython_files.tar"

clear_cache=0
home_dir=$(pwd)

if [ -f "$CACHE_File" ] && [ "$USE_CACHE" ] && [ -d "$PYX_CACHE_DIR" ]; then

    echo "Cache available - checking pyx diff"

    for i in ${pyx_file_list}
    do
            diff=`diff -u $i $PYX_CACHE_DIR${i}`
            if [[ $? -ne 0 ]]
            then
                    echo "${i##*/} can't be diffed; probably not in cache"
                    clear_cache=1
            fi
            if [[ ! -z $diff ]]
            then
                    echo "${i##*/} has changed:"
                    echo $diff
                    clear_cache=1
            fi
    done

    if [ "$TRAVIS_PULL_REQUEST" == "false" ]
    then
        echo "Not a PR"
        # Uncomment next 2 lines to turn off cython caching not in a PR
        # echo "Non PR cython caching is disabled"
        # clear_cache=1
    else
        echo "In a PR"
        # Uncomment next 2 lines to turn off cython caching in a PR
        # echo "PR cython caching is disabled"
        # clear_cache=1
    fi

fi

if [ $clear_cache -eq 1 ] && [ "$USE_CACHE" ]
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
