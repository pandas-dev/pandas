#!/bin/bash

# This must be sourced by .travis.yml, so any envars exported here will
# be available to the rest of the build stages

# - computes a hash based on the cython files in the codebade
# - retrieves the decrypted key if any for all whitelisted forks
# - checks whether the user optd int to use the cache
# - if so, check for availablity of cache files on the server, based on hash
# - set envars to control what the following scripts do

# at most one of these will decrypt, so the end result is that $STORE_KEY
# either holds a single key or does not
export STORE_KEY="$STORE_KEY0""$STORE_KEY1""$STORE_KEY2""$STORE_KEY3""$STORE_KEY4"
export STORE_KEY="$STORE_KEY""$STORE_KEY5""$STORE_KEY6""$STORE_KEY7"

export CYTHON_HASH=$(find pandas | grep -P '\.(pyx|pxd)$'  | sort \
     | while read N; do echo $(tail -n+1 $N | md5sum ) ;done | md5sum| cut -d ' ' -f 1)

export CYTHON_HASH=$CYTHON_HASH-$TRAVIS_PYTHON_VERSION

# where the cache files live on the server
export CACHE_FILE_URL="https://cache27-pypandas.rhcloud.com/static/$STORE_KEY/$CYTHON_HASH.zip"
export VENV_FILE_URL="https://cache27-pypandas.rhcloud.com/static/$STORE_KEY/venv-$TRAVIS_PYTHON_VERSION.zip"
export CACHE_FILE_STORE_URL="https://cache27-pypandas.rhcloud.com/store/$STORE_KEY"

echo "Hashing:"
find pandas | grep -P '\.(pyx|pxd)$'
echo "Key: $CYTHON_HASH"

export CACHE_FILE_AVAILABLE=false
export VENV_FILE_AVAILABLE=false
export PLEASE_TRAVIS_FASTER=false

# check whether the user opted in to use the cache via commit message
if [ x"$(git log --format='%s' -n 1 | grep PLEASE_TRAVIS_FASTER | wc -l)" != x"0" ]; then
    export PLEASE_TRAVIS_FASTER=true
fi;
if [ x"$(git log --format='%s' -n 1 | grep PTF | wc -l)" != x"0" ]; then
    export PLEASE_TRAVIS_FASTER=true
fi;

if $PLEASE_TRAVIS_FASTER; then

    # check whether the files exists on the server
    curl -s -f -I "$CACHE_FILE_URL"  #  silent, don;t expose key
    if [ x"$?" == x"0" ] ; then
        export CACHE_FILE_AVAILABLE=true;
    fi


    curl -s -f -I "$VENV_FILE_URL" #  silent, don;t expose key
    if [ x"$?" == x"0" ] ; then
        export VENV_FILE_AVAILABLE=true;
    fi

    # the pandas build cache machinery needs this set, and the directory created
    export BUILD_CACHE_DIR="/tmp/build_cache"
    mkdir "$BUILD_CACHE_DIR"
fi;

# debug
echo "PLEASE_TRAVIS_FASTER=$PLEASE_TRAVIS_FASTER"
echo "CACHE_FILE_AVAILABLE=$CACHE_FILE_AVAILABLE"
echo "VENV_FILE_AVAILABLE=$VENV_FILE_AVAILABLE"

true
