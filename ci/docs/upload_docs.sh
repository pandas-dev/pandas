#!/bin/bash
set -e

echo "inside $0"

if [ "${DOC}" ] && [ "${TRAVIS_PULL_REQUEST}" = "false" ]; then

    echo ########################
    echo # Create and send docs #
    echo ########################

    cd /tmp/doc/build/html
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

    git remote remove origin
    git remote add origin "https://${PANDAS_GH_TOKEN}@github.com/pandas-dev/pandas-docs-travis.git"
    git fetch origin
    git remote -v

    git push origin gh-pages -f
else
    echo "[skipping doc upload]"
fi

exit 0
