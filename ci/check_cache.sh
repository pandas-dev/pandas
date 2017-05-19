#!/bin/bash

# currently not used
# script to make sure that cache is clean
# Travis CI now handles this

if [ "$TRAVIS_PULL_REQUEST" == "false" ]
then
    echo "Not a PR: checking for changes in ci/ from last 2 commits"
    git diff HEAD~2 --numstat | grep -E "ci/"
    ci_changes=$(git diff HEAD~2 --numstat | grep -E "ci/"| wc -l)
else
    echo "PR: checking for changes in ci/ from last 2 commits"
    git fetch origin pull/${TRAVIS_PULL_REQUEST}/head:PR_HEAD
    git diff PR_HEAD~2 --numstat | grep -E "ci/"
    ci_changes=$(git diff PR_HEAD~2 --numstat | grep -E "ci/"| wc -l)
fi

CACHE_DIR="$HOME/.cache/"
CCACHE_DIR="$HOME/.ccache/"

if [ $ci_changes -ne 0 ]
then
    echo "Files have changed in ci/ deleting all caches"
    rm -rf "$CACHE_DIR"
    rm -rf "$CCACHE_DIR"
fi
