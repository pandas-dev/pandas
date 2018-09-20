#!/bin/bash
set -e

echo "inside $0"

if [ "$DOC" ]; then
    cd /tmp/doc

    echo "[linting docs]"

    echo './make.py lint_log --log-file=doc-build.log'
    ./make.py lint_log --log-file=doc-build.log
else
    echo "[skipping doc lint]"
fi
