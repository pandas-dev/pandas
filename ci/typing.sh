#!/bin/bash

echo "inside $0"

source activate pandas

RET=0

if [ "$TYPING" ]; then

    echo "Typing  *.py"
    mypy \
        pandas/core/base.py
    if [ $? -ne "0" ]; then
        RET=1
    fi
    echo "Typing *.py DONE"

else
    echo "NOT checking typing"
fi

exit $RET
