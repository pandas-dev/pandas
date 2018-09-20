#!/bin/bash
set -e

cd "$TRAVIS_BUILD_DIR"
echo "inside $0"

if [ "$DOC" ]; then

    echo "[building docs]"

    source activate pandas

    mv "$TRAVIS_BUILD_DIR"/doc /tmp
    mv "$TRAVIS_BUILD_DIR/LICENSE" /tmp  # included in the docs.
    cd /tmp/doc

    echo './make.py 2>&1 | tee doc-build.log'
    ./make.py 2>&1 | tee doc-build.log
else
   echo "[skipping docs]"
fi

exit 0
