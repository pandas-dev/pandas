#!/bin/bash

echo "inside $0"

if [ -n "$LOCALE_OVERRIDE" ]; then
    export LC_ALL="$LOCALE_OVERRIDE";
    echo "Setting LC_ALL to $LOCALE_OVERRIDE"
    curdir="$(pwd)"
    cd /tmp
    pycmd='import pandas; print("pandas detected console encoding: %s" % pandas.get_option("display.encoding"))'
    python -c "$pycmd"
    cd "$curdir"
fi

echo nosetests -v --exe -w /tmp -A "$NOSE_ARGS" pandas --with-xunit --xunit-file=/tmp/nosetests.xml
nosetests -v --exe -w /tmp -A "$NOSE_ARGS" pandas --with-xunit --xunit-file=/tmp/nosetests.xml
