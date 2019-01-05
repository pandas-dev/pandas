#!/bin/bash

set -e

if [ "${TRAVIS_OS_NAME}" != "linux" ]; then
   echo "not doing build_docs on non-linux"
   exit 0
fi

cd "$TRAVIS_BUILD_DIR"/doc
echo "inside $0"

if [ "$DOC" ]; then

    echo "Will build docs"

    echo ###############################
    echo # Log file for the doc build  #
    echo ###############################

    echo ./make.py
    ./make.py

    echo ########################
    echo # Create and send docs #
    echo ########################

    echo "Only uploading docs when TRAVIS_PULL_REQUEST is 'false'"
    echo "TRAVIS_PULL_REQUEST: ${TRAVIS_PULL_REQUEST}"

    if [ "${TRAVIS_PULL_REQUEST}" == "false" ]; then
        cd build/html
        git config --global user.email "pandas-docs-bot@localhost.foo"
        git config --global user.name "pandas-docs-bot"

        # create the repo
        git init

        touch README
        git add README
        git commit -m "Initial commit" --allow-empty
        git branch gh-pages
        git checkout gh-pages
        touch .nojekyll
        git add --all .
        git commit -m "Version" --allow-empty

        git remote add origin "https://${PANDAS_GH_TOKEN}@github.com/pandas-dev/pandas-docs-travis.git"
        git fetch origin
        git remote -v

        git push origin gh-pages -f
    fi
fi

exit 0
