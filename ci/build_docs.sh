#!/bin/bash


cd "$TRAVIS_BUILD_DIR"

git show --pretty="format:" --name-only HEAD~5.. --first-parent | grep -P "rst|txt|doc"

if [ "$?" != "0" ]; then
    echo "Skipping doc build, none were modified"
    # nope, skip docs build
    exit 0
fi


if [ x"$DOC_BUILD" != x"" ]; then
    # we're running network tests, let's build the docs in the meantim
    echo "Will build docs"
    pip install sphinx==1.1.3 ipython==1.1.0

    mv "$TRAVIS_BUILD_DIR"/doc /tmp
    cd /tmp/doc

    rm /tmp/doc/source/api.rst # no R
    rm /tmp/doc/source/r_interface.rst # no R

    echo ############################### > /tmp/doc.log
    echo # Log file for the doc build  # > /tmp/doc.log
    echo ############################### > /tmp/doc.log
    echo "" > /tmp/doc.log
    echo -e "y\n" | ./make.py --no-api 2>&1

    cd /tmp/doc/build/html
    git config --global user.email "pandas-docs-bot@localhost.foo"
    git config --global user.name "pandas-docs-bot"

    git init
    touch README
    git add README
    git commit -m "Initial commit" --allow-empty
    git branch gh-pages
    git checkout gh-pages
    touch .nojekyll
    git add --all .
    git commit -m "Version" --allow-empty
    git remote add origin https://$GH_TOKEN@github.com/pandas-docs/pandas-docs-travis
    git push origin gh-pages -f
fi

exit 0
