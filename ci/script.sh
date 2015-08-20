#!/bin/bash

echo "inside $0"

source activate pandas

if [ -n "$LOCALE_OVERRIDE" ]; then
    export LC_ALL="$LOCALE_OVERRIDE";
    echo "Setting LC_ALL to $LOCALE_OVERRIDE"

    pycmd='import pandas; print("pandas detected console encoding: %s" % pandas.get_option("display.encoding"))'
    python -c "$pycmd"
fi

if [ "$BUILD_TEST" ]; then
    echo "We are not running nosetests as this is simply a build test."
else
    echo nosetests --exe -A "$NOSE_ARGS" pandas --doctest-tests --with-xunit --xunit-file=/tmp/nosetests.xml
    nosetests --exe -A "$NOSE_ARGS" pandas --doctest-tests --with-xunit --xunit-file=/tmp/nosetests.xml
fi

RET="$?"

exit "$RET"
