#!/bin/bash

echo "inside $0"

source activate pandas

if [ -n "$LOCALE_OVERRIDE" ]; then
    export LC_ALL="$LOCALE_OVERRIDE";
    echo "Setting LC_ALL to $LOCALE_OVERRIDE"

    pycmd='import pandas; print("pandas detected console encoding: %s" % pandas.get_option("display.encoding"))'
    python -c "$pycmd"
fi

# conditionally build and upload docs to GH/pandas-docs/pandas-docs/travis
"$TRAVIS_BUILD_DIR"/ci/build_docs.sh 2>&1 > /tmp/doc.log &
# doc build log will be shown after tests

if [ "$BUILD_TEST" ]; then
    echo "We are not running nosetests as this is simply a build test."
else
    echo nosetests --exe -A "$NOSE_ARGS" pandas --with-xunit --xunit-file=/tmp/nosetests.xml
    nosetests --exe -A "$NOSE_ARGS" pandas --with-xunit --xunit-file=/tmp/nosetests.xml
fi

if [ "$INSTALL_TEST" ]; then
    echo "Starting installation test."
    conda uninstall cython || exit 1
    python "$TRAVIS_BUILD_DIR"/setup.py sdist --formats=zip,gztar || exit 1
    pip install "$TRAVIS_BUILD_DIR"/dist/*tar.gz || exit 1
    nosetests --exe -A "$NOSE_ARGS" pandas/tests/test_series.py --with-xunit --xunit-file=/tmp/nosetests_install.xml
fi
RET="$?"

# wait until subprocesses finish (build_docs.sh)
wait

exit "$RET"
