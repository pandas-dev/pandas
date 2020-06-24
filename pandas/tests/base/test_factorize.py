import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize("sort", [True, False])
def test_factorize(index_or_series_obj, sort):
    obj = index_or_series_obj
    result_codes, result_uniques = obj.factorize(sort=sort)

    constructor = pd.Index
    if isinstance(obj, pd.MultiIndex):
        constructor = pd.MultiIndex.from_tuples
    expected_uniques = constructor(obj.unique())

    if sort:
        expected_uniques = expected_uniques.sort_values()

    # construct an integer ndarray so that
    # `expected_uniques.take(expected_codes)` is equal to `obj`
    expected_uniques_list = list(expected_uniques)
    expected_codes = [expected_uniques_list.index(val) for val in obj]
    expected_codes = np.asarray(expected_codes, dtype=np.intp)

    tm.assert_numpy_array_equal(result_codes, expected_codes)
    tm.assert_index_equal(result_uniques, expected_uniques)
