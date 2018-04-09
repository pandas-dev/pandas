#!/bin/bash

echo "inside $0"


source activate pandas
cd "$TRAVIS_BUILD_DIR"

RET=0

if [ "$DOCTEST" ]; then

    echo "Running doctests"

    # top-level reshaping functions
    pytest --doctest-modules -v \
        pandas/core/reshape/concat.py \
        pandas/core/reshape/pivot.py \
        pandas/core/reshape/reshape.py \
        pandas/core/reshape/tile.py

    if [ $? -ne "0" ]; then
        RET=1
    fi

    # top-level io functionality
    pytest --doctest-modules -v pandas/io/* -k "read_"

    if [ $? -ne "0" ]; then
        RET=1
    fi

    # DataFrame docstrings
    pytest --doctest-modules -v \
        pandas/core/frame.py -k"-agg -apply -applymap -assign -eval -isin -itertuples -join -merge -nlargest -nsmallest -nunique -pivot -pivot_table -quantile -query -reindex -reindex_axis -replace -round -select_dtypes -set_index -stack -to_dict -to_excel -to_stata"

    if [ $? -ne "0" ]; then
        RET=1
    fi

    # Series docstrings
    pytest --doctest-modules -v \
        pandas/core/series.py -k"-agg -apply -map -nlargest -nonzero -nsmallest -reindex -replace -searchsorted -to_dict -to_excel"

    if [ $? -ne "0" ]; then
        RET=1
    fi

else
    echo "NOT running doctests"
fi

exit $RET
