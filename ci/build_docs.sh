#!/bin/bash

if [ "${TRAVIS_OS_NAME}" != "linux" ]; then
   echo "not doing build_docs on non-linux"
   exit 0
fi

cd "$TRAVIS_BUILD_DIR"
echo "inside $0"

git show --pretty="format:" --name-only HEAD~5.. --first-parent | grep -P "rst|txt|doc"

if [ "$?" != "0" ]; then
    echo "Skipping doc build, none were modified"
    # nope, skip docs build
    exit 0
fi


if [ x"$DOC_BUILD" != x"" ]; then

    echo "Will build docs"

    source activate pandas
    conda install -n pandas -c r r rpy2 --yes

    time sudo apt-get $APT_ARGS install dvipng texlive-latex-base texlive-latex-extra

    mv "$TRAVIS_BUILD_DIR"/doc /tmp
    cd /tmp/doc

    echo ###############################
    echo # Log file for the doc build  #
    echo ###############################

    echo ./make.py
    ./make.py

    echo ########################
    echo # Create and send docs #
    echo ########################

    cd /tmp/doc/build/html
    git config --global user.email "pandas-docs-bot@localhost.foo"
    git config --global user.name "pandas-docs-bot"
    git config --global credential.helper cache

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
    git remote add origin "https://$GH_TOKEN@github.com/pandas-docs/pandas-docs-travis"
    git push origin gh-pages -f
fi

exit 0
