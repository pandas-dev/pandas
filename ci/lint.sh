#!/bin/bash

echo "inside $0"

source activate pandas

RET=0

if [ "$LINT" ]; then
    echo "Linting"
    for path in 'core' 'indexes' 'types' 'formats' 'io' 'stats' 'compat' 'sparse' 'tools' 'tseries' 'tests' 'computation' 'util'
    do
        echo "linting -> pandas/$path"
        flake8 pandas/$path --filename '*.py'
        if [ $? -ne "0" ]; then
            RET=1
        fi
    done
else
    echo "NOT Linting"
fi

exit $RET
