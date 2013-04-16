#!/bin/bash

# If envars.sh determined we're running  in an authorized fork
# and the user opted in to the network cache,and that cached versions
# are available on the cache server, download and deploy the cached
# files to the local filesystem

echo "inside $0"

# overview
sudo apt-get update $APT_ARGS # run apt-get update for all versions

if  $PLEASE_TRAVIS_FASTER ; then
    echo "Faster? well... I'll try."

    if $CACHE_FILE_AVAILABLE ;  then
        echo retrieving "$CACHE_FILE_URL";

        wget -q "$CACHE_FILE_URL" -O "/tmp/_$CYTHON_HASH.zip";
        unzip $ZIP_FLAGS /tmp/_"$CYTHON_HASH.zip" -d "$BUILD_CACHE_DIR";
        rm -f /tmp/_"$CYTHON_HASH.zip"
        # copy cythonized c files over
        cp -R "$BUILD_CACHE_DIR"/pandas/* pandas/
        # mkdir build/
        # ls -l build
        # mkdir build/temp.linux-x86-$TRAVIS_PYTHON_VERSION
        # mkdir build/temp.linux-x86-$TRAVIS_PYTHON_VERSION/pandas/
        # touch build/temp.linux-x86-$TRAVIS_PYTHON_VERSION/pandas/lib.o

    fi;
    echo "VENV_FILE_AVAILABLE=$VENV_FILE_AVAILABLE"
    if $VENV_FILE_AVAILABLE ; then
        echo "getting venv"
        wget -q $VENV_FILE_URL -O  "/tmp/venv.zip";
        sudo unzip $ZIP_FLAGS -o /tmp/venv.zip -d "/";
        sudo chown travis -R "$VIRTUAL_ENV"
        rm -f /tmp/_"$CYTHON_HASH.zip"
    fi;
fi

true # never fail because bad things happened here
