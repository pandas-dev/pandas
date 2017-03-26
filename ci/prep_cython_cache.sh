#!/bin/bash

ls "$HOME/.cache/"

PYX_CACHE_DIR="$HOME/.cache/pyxfiles"
pyx_file_list=`find ${TRAVIS_BUILD_DIR} -name "*.pyx" -o -name "*.pxd" -o -name "*.pxi.in"`
pyx_cache_file_list=`find ${PYX_CACHE_DIR} -name "*.pyx" -o -name "*.pxd" -o -name "*.pxi.in"`

CACHE_File="$HOME/.cache/cython_files.tar"

# Clear the cython cache 0 = NO, 1 = YES
clear_cache=0

pyx_files=`echo "$pyx_file_list" | wc -l`
pyx_cache_files=`echo "$pyx_cache_file_list" | wc -l`

if [[ pyx_files -ne pyx_cache_files ]]
then
        echo "Different number of pyx files"
        clear_cache=1
fi

home_dir=$(pwd)

if [ -f "$CACHE_File" ] && [ -z "$NOCACHE" ] && [ -d "$PYX_CACHE_DIR" ]; then

    echo "Cache available - checking pyx diff"

    for i in ${pyx_file_list}
    do
            diff=`diff -u $i $PYX_CACHE_DIR${i}`
            if [[ $? -eq 2 ]]
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

if [ $clear_cache -eq 0 ] && [ -z "$NOCACHE" ]
then
    # No and nocache is not set
    echo "Will reuse cached cython file"
    cd /
    tar xvmf $CACHE_File
    cd $home_dir
else
    echo "Rebuilding cythonized files"
    echo "No cache = $NOCACHE"
    echo "Clear cache (1=YES) = $clear_cache"
fi


exit 0
