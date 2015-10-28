#!/bin/bash

# build the distribution
LAST=`git tag --sort version:refname | tail -1`

echo "Building distribution for: $LAST"
git checkout $LAST

read -p "Ok to continue (y/n)? " answer
case ${answer:0:1} in
    y|Y )
        echo "Building distribution"
        python setup.py clean
        python setup.py build_ext --inplace
        python setup.py sdist --formats=zip,gztar
    ;;
    * )
        echo "Not building distribution"
    ;;
esac
