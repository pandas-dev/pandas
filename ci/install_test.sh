#!/bin/bash

echo "inside $0"

if [ "$INSTALL_TEST" ]; then
    source activate pandas
    echo "Starting installation test."
    conda uninstall cython || exit 1
    python "$TRAVIS_BUILD_DIR"/setup.py sdist --formats=zip,gztar || exit 1
    pip install "$TRAVIS_BUILD_DIR"/dist/*tar.gz || exit 1
    nosetests --exe -A "$NOSE_ARGS" pandas/tests/test_series.py --with-xunit --xunit-file=/tmp/nosetests_install.xml
else
    echo "Skipping installation test."
fi
RET="$?"

exit "$RET"
