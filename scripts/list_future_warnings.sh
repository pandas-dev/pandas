#!/bin/bash

# Check all future warnings in Python files, and report them with the version
# where the FutureWarning was added.
#
# This is useful to detect features that have been deprecated, and should be
# removed from the code. For example, if a line of code contains:
#
#     warning.warn('Method deprecated', FutureWarning, stacklevel=2)
#
# Which is released in Pandas 0.20.0, then it is expected that the method
# is removed before releasing Pandas 0.24.0, including the warning. If it
# is not, this script will list this line, with the version 0.20.0, which
# will make it easy to detect that it had to be removed.
#
# In some cases this script can return false positives, for example in files
# where FutureWarning is used to detect deprecations, or similar. The EXCLUDE
# variable can be used to ignore files that use FutureWarning, but do not
# deprecate functionality.
#
# Usage:
#
#     $ ./list_future_warnings.sh

EXCLUDE="^pandas/tests/|"  # tests validate that FutureWarnings are raised
EXCLUDE+="^pandas/util/_decorators.py$|"  # generic deprecate function that raises warning
EXCLUDE+="^pandas/util/_depr_module.py$|"  # generic deprecate module that raises warnings
EXCLUDE+="^pandas._testing.py$|" # contains function to evaluate if warning is raised
EXCLUDE+="^pandas/io/parsers.py$"  # implements generic deprecation system in io reading

BASE_DIR="$(dirname $0)/.."
cd $BASE_DIR
FILES=`grep -RIl "FutureWarning" pandas/* | grep -vE "$EXCLUDE"`
OUTPUT=()
IFS=$'\n'

for FILE in $FILES; do
    FILE_LINES=`git blame -sf $FILE | grep FutureWarning | tr -s " " | cut -d " " -f1,3`
    for FILE_LINE in $FILE_LINES; do
        TAG=$(git tag --contains $(echo $FILE_LINE | cut -d" " -f1) | head -n1)
        OUTPUT_ROW=`printf "%-14s %-16s %s" ${TAG:-"(not released)"} $FILE_LINE $FILE`
        OUTPUT+=($OUTPUT_ROW)
    done
done

printf "%s\n" "${OUTPUT[@]}" | sort -V
