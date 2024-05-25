import pytest
from pandas.core.array_algos.find_min_subset_sum import find_min_subset_sum

test_cases = [
    # Exact match no restrictions
    ([5, 1, 8, 3, 12, 9], 15, False,
     {'subset': [3, 12], 'sum': 15}),

    # Exact match no restrictions with duplicated values
    ([5, 2, 0, 3, 5], 15, False,
     {'subset': [5, 2, 3, 5], 'sum': 15}),

    # Closest match constraint
    ([5, 8, 3, 12, 9, 1], 16, True,
     {'subset': [5, 8, 3], 'sum': 16}),

    # Closest match onstraint with duplicated values and zero
    ([5, 2, 0, 3, 5], 16, True,
     {'subset': [5, 2, 3, 5], 'sum': 15}),

    # No exact match, return empty array
    ([5, 1, 8, 3, 12, 9], 100, False,
     {'subset': [], 'sum': 0}),

    # Empty array input
    ([], 10, False,
     {'subset': [], 'sum': 0}),

    # Single element match
    ([10], 10, False,
     {'subset': [10], 'sum': 10}),

    # Single element less than targets closest match
    ([5], 10, True,
     {'subset': [5], 'sum': 5}),
]


@pytest.mark.parametrize('arr, target, closest, expected', test_cases)
def test_find_min_subset_sum(arr, target, closest, expected):
    assert find_min_subset_sum(arr, target, closest) == expected


# Test cases for ValueError scenarios
invalid_arr_cases = [
    ([1, 2, 'a'], 10),  # Non-integer value in arr
    (123, 10),           # Non-list value for arr
    ([-1, 2, 3], 10),    # Negative integer in arr
]

invalid_target_cases = [
    ([1, 2, 3], 'abc'),  # Non-integer value for target
]

@pytest.mark.parametrize('arr, target', invalid_arr_cases + invalid_target_cases)
def test_invalid_inputs(arr, target):
    with pytest.raises(ValueError):
        find_min_subset_sum(arr, target)
