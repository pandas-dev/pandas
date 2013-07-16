#!/bin/bash

echo "inside $0"

if [ x"$LOCALE_OVERRIDE" != x"" ]; then
    export LC_ALL="$LOCALE_OVERRIDE";
    echo "Setting LC_ALL to $LOCALE_OVERRIDE"
    (cd /; python -c 'import pandas; print("pandas detected console encoding: %s" % pandas.get_option("display.encoding"))')

fi

echo nosetests --exe -w /tmp -A "$NOSE_ARGS" pandas --with-xunit --xunit-file=/tmp/nosetests.xml
nosetests --exe -w /tmp -A "$NOSE_ARGS" pandas --with-xunit --xunit-file=/tmp/nosetests.xml
