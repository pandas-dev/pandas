#!/bin/bash

echo "inside $0"

if [ -n "$LOCALE_OVERRIDE" ]; then
    export LC_ALL="$LOCALE_OVERRIDE";
    echo "Setting LC_ALL to $LOCALE_OVERRIDE"
    (cd /; python -c 'import pandas; print("pandas detected console encoding: %s" % pandas.get_option("display.encoding"))')

fi

echo nosetests --exe -w /tmp -A "$NOSE_ARGS" pandas --show-skipped
nosetests --exe -w /tmp -A "$NOSE_ARGS" pandas --show-skipped
