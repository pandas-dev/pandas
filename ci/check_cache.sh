#!/bin/bash

git diff HEAD~2 --numstat | grep -E "ci/"
ci_changes=$(git diff HEAD~2 --numstat | grep -E "ci/"| wc -l)

MINICONDA_DIR="$HOME/miniconda/"
CACHE_DIR="$HOME/.cache/"
CCACHE_DIR="$HOME/.ccache/"

if [ $ci_changes -ne 0 ]
then
    echo "CI has been changed in the last two commits deleting all caches"
    rm -rf "$MINICONDA_DIR"
    rm -rf "$CACHE_DIR"
    rm -rf "$CCACHE_DIR"
fi