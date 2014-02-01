#!/bin/bash

if [ "$IRON_TOKEN" ]; then
    sudo apt-get $APT_ARGS install ccache p7zip-full
    # iron_cache, pending py3 fixes upstream
    pip install -I --allow-external --allow-insecure git+https://github.com/y-p/iron_cache_python.git


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
fi

exit 0
