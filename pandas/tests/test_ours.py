import sys
import os
current_script_dir = os.path.dirname(os.path.abspath(__file__))
local_pandas_path = os.path.abspath(os.path.join(current_script_dir, "../../"))
sys.path.insert(0, local_pandas_path)

import numpy as np
import pytest
from pandas.core.nanops import get_corr_func
from pandas.core.nanops import branch_coverage_get_corr_func

def print_branch_coverage(branch_coverage):
    trues = 0
    for key in branch_coverage:
        print(f"{key}: {branch_coverage[key]}")
        if branch_coverage[key] is True:
            trues += 1
    
    total = len(branch_coverage)
    coverage = (trues / total) * 100
    print(f"Branch coverage: {coverage}%")

def test_get_corr_func():
    # Test data
    a = np.array([1, 2, 3, 4, 5])
    b = np.array([5, 4, 3, 2, 1])

    def test_kendall():
        corr_func = get_corr_func("kendall")
        assert isinstance(corr_func(a, b), float)

    def test_spearman():
        corr_func = get_corr_func("spearman")
        assert isinstance(corr_func(a, b), float)

    def test_pearson():
        corr_func = get_corr_func("pearson")
        assert isinstance(corr_func(a, b), float)

    def test_callable_method():
        def custom_corr(a, b):
            return 0.5
        corr_func = get_corr_func(custom_corr)
        assert corr_func(a, b) == 0.5

    def test_unknown_method():
        with pytest.raises(ValueError):
            get_corr_func("unknown")

    test_kendall()
    test_spearman()
    test_pearson()
    test_callable_method()
    test_unknown_method()

# test_get_corr_func()
# print_branch_coverage(branch_coverage_get_corr_func)


from pandas.core.indexes.range import RangeIndex
from pandas.core.indexes.range import branch_coverage_sort_values

def test_sort_values():
    def test_return_indexer_false_ascending():
        ri = RangeIndex(0, 10, 1)
        result, indexer = ri.sort_values(return_indexer=True, ascending=True)
        assert indexer is not None, "Indexer should not be None when return_indexer is True"

    def test_return_indexer_true_ascending():
        ri = RangeIndex(0, 10, -1)
        result, indexer = ri.sort_values(return_indexer=True, ascending=True)
        assert indexer is not None, "Indexer should not be None when return_indexer is True"

    def test_return_indexer_true_descending():
        ri = RangeIndex(0, 10, 1)
        result, indexer = ri.sort_values(return_indexer=True, ascending=False)
        assert indexer is not None, "Indexer should not be None when return_indexer is True"

    def test_return_indexer_false_descending():
        ri = RangeIndex(0, 10, -1)
        result, indexer = ri.sort_values(return_indexer=True, ascending=False)
        assert indexer is not None, "Indexer should not be None when return_indexer is True"


    def test_sort_values_with_key():
        ri = RangeIndex(-5, 5, 1)
        # Using a key function that sorts based on the absolute difference from 5
        result, indexer = ri.sort_values(return_indexer=True, ascending=True, key=lambda x: abs(x - 5))
        expected_result = [5, 4, 6, 3, 7, 2, 8, 1, 9, 0, -5]
        assert list(result) == expected_result, "Result should be sorted based on the absolute difference from 5"
        assert indexer is not None, "Indexer should not be None when return_indexer is True"

    # test_sort_values_with_key()
    test_return_indexer_false_ascending()
    test_return_indexer_true_ascending()
    test_return_indexer_true_descending()
    test_return_indexer_false_descending()

test_sort_values()
print_branch_coverage(branch_coverage_sort_values)
