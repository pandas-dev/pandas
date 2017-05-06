#!/bin/bash

# build the distribution
LAST=`git tag --sort version:refname | grep -v rc | tail -1`

echo "Building distribution for: $LAST"
git checkout $LAST

read -p "Ok to continue (y/n)? " answer
case ${answer:0:1} in
    y|Y )
        echo "Building distribution"
        ./build_dist_for_release.sh
    ;;
    * )
        echo "Not building distribution"
    ;;
esac
