#!/bin/bash

# build the distribution
LAST=`git tag --sort version:refname | grep -v rc | tail -1`

echo "Building distribution for: $LAST"
git checkout $LAST

read -p "Ok to continue (y/n)? " answer
case ${answer:0:1} in
    y|Y )
        echo "Building distribution"
        rm -rf dist
        git clean -xfd
        python setup.py clean
        python setup.py cython
        python setup.py sdist --formats=gztar
    ;;
    * )
        echo "Not building distribution"
    ;;
esac
