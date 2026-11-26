import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize(
    "values,dtype",
    [
        ([2, 1, 2, None, None], "Int32"),
        ([2**63 + 1, 2**63, 2**63 + 1, None, None], "UInt64"),
        ([2.5, 1.5, 2.5, None, None], "Float32"),
        ([True, False, True, None, None], "boolean"),
        (["b", "a", "b", None, None], "string"),
        (
            ["b", "a", "b", None, None],
            pd.CategoricalDtype(["b", "a", "unused"], ordered=True),
        ),
    ],
)
@pytest.mark.parametrize("nlevels", [1, 2])
@pytest.mark.parametrize("sort", [True, False])
@pytest.mark.parametrize("use_na_sentinel", [True, False])
def test_factorize_preserves_level_dtypes(
    values, dtype, nlevels, sort, use_na_sentinel
):
    # GH#62337
    arrays = [pd.array(values, dtype=dtype)]
    if nlevels == 2:
        arrays.append(pd.array([1, 2, 1, 0, 0], dtype="Int16"))
    index = pd.MultiIndex.from_arrays(arrays, names=list("ab")[:nlevels])

    codes, uniques = index.factorize(sort=sort, use_na_sentinel=use_na_sentinel)

    expected_codes = np.array([0, 1, 0, 2, 2], dtype=np.intp)
    positions = [0, 1, 3]
    if sort:
        expected_codes = np.array([1, 0, 1, 2, 2], dtype=np.intp)
        positions = [1, 0, 3]
    expected = index.take(positions)

    tm.assert_numpy_array_equal(codes, expected_codes)
    tm.assert_index_equal(uniques, expected, exact=True)
    tm.assert_index_equal(uniques.take(codes), index, exact=True)


@pytest.mark.parametrize("values", [[], [None, None]])
@pytest.mark.parametrize("sort", [True, False])
def test_factorize_empty_or_all_na_preserves_levels(values, sort):
    # GH#62337
    index = pd.MultiIndex.from_arrays([pd.array(values, dtype="Int32")], names=["col"])

    codes, uniques = pd.factorize(index, sort=sort)

    expected_codes = np.zeros(len(index), dtype=np.intp)
    expected = index[:1]
    tm.assert_numpy_array_equal(codes, expected_codes)
    tm.assert_index_equal(uniques, expected, exact=True)


@pytest.mark.parametrize("sort", [True, False])
def test_factorize_many_levels(sort):
    # GH#62337: the cartesian product of the levels does not fit in int64.
    index = pd.MultiIndex.from_arrays([pd.array([2, 1, 2, None], dtype="Int32")] * 65)

    codes, uniques = index.factorize(sort=sort)

    expected_codes = np.array([1, 0, 1, 2] if sort else [0, 1, 0, 2], dtype=np.intp)
    expected = index.take([1, 0, 3] if sort else [0, 1, 3])
    tm.assert_numpy_array_equal(codes, expected_codes)
    tm.assert_index_equal(uniques, expected, exact=True)
