"""
these are systematically testing all of the args to value_counts
with different size combinations. This is to ensure stability of the sorting
and proper parameter handling
"""

from itertools import product

import numpy as np
import pytest

from pandas import DataFrame, Grouper, MultiIndex, Series, date_range, to_datetime
import pandas._testing as tm


# our starting frame
def seed_df(seed_nans, n, m):
    np.random.seed(1234)
    days = date_range("2015-08-24", periods=10)

    frame = DataFrame(
        {
            "1st": np.random.choice(list("abcd"), n),
            "2nd": np.random.choice(days, n),
            "3rd": np.random.randint(1, m + 1, n),
        }
    )

    if seed_nans:
        frame.loc[1::11, "1st"] = np.nan
        frame.loc[3::17, "2nd"] = np.nan
        frame.loc[7::19, "3rd"] = np.nan
        frame.loc[8::19, "3rd"] = np.nan
        frame.loc[9::19, "3rd"] = np.nan

    return frame


# create input df, keys, and the bins
binned = []
ids = []
for seed_nans in [True, False]:
    for n, m in product((100, 1000), (5, 20)):
        df = seed_df(seed_nans, n, m)
        bins = None, np.arange(0, max(5, df["3rd"].max()) + 1, 2)
        keys = "1st", "2nd", ["1st", "2nd"]
        for k, b in product(keys, bins):
            binned.append((df, k, b, n, m))
            ids.append(f"{k}-{n}-{m}-{seed_nans} ")


@pytest.mark.slow
@pytest.mark.parametrize("df, keys, bins, n, m", binned, ids=ids)
@pytest.mark.parametrize("isort", [True, False])
@pytest.mark.parametrize("normalize", [True, False])
@pytest.mark.parametrize("sort", [True, False])
@pytest.mark.parametrize("ascending", [True, False])
@pytest.mark.parametrize("dropna", [True, False])
def test_series_groupby_value_counts(
    df, keys, bins, n, m, isort, normalize, sort, ascending, dropna
):
    def rebuild_index(df):
        arr = list(map(df.index.get_level_values, range(df.index.nlevels)))
        df.index = MultiIndex.from_arrays(arr, names=df.index.names)
        return df

    kwargs = dict(
        normalize=normalize, sort=sort, ascending=ascending, dropna=dropna, bins=bins
    )

    print(f"{df=}")
    gr = df.groupby(keys, sort=isort)
    left = gr["3rd"].value_counts(**kwargs)
    left.index.names = left.index.names[:-1] + ["3rd"]

    # gr = df.groupby(keys, sort=isort)
    right = gr["3rd"].apply(Series.value_counts, **kwargs)
    right.index.names = right.index.names[:-1] + ["3rd"]
    print(f"{left=}")
    print(f"{right=}")

    # have to sort on index because of unstable sort on values
    left, right = map(rebuild_index, (left, right))  # xref GH9212
    # have to ignore 0 counts to be consistent with individual column value_counts
    left = left[left.astype(bool)]
    right = right[right.astype(bool)]
    tm.assert_series_equal(left.sort_index(), right.sort_index())


def test_groubpy_value_counts_bins():
    # GH32471
    BINS = [0, 20, 80, 100]
    values = [
        [0, 5, 0],
        [1, 5, 100],
        [0, 5, 100],
        [2, 5, 0],
        [3, 6, 100],
        [3, 5, 100],
        [1, 5, 100],
    ]
    df = DataFrame(values, columns=["key1", "key2", "score"])
    result = df.groupby(["key1", "key2"])["score"].value_counts(bins=BINS)
    print(f"{result=}")
    print(
        df.groupby(["key1", "key2"])["score"].apply(
            Series.value_counts,
            bins=BINS,
            sort=True,
            normalize=True,
            ascending=True,
            dropna=True,
        )
    )

    result.sort_index(inplace=True)
    expected = Series(
        [1, 0, 1, 0, 0, 2, 1, 0, 0, 0, 0, 1, 0, 0, 1], result.index, name="score"
    )
    tm.assert_series_equal(result, expected)


def test_series_groupby_value_counts_with_grouper():
    # GH28479
    df = DataFrame(
        {
            "Timestamp": [
                1565083561,
                1565083561 + 86400,
                1565083561 + 86500,
                1565083561 + 86400 * 2,
                1565083561 + 86400 * 3,
                1565083561 + 86500 * 3,
                1565083561 + 86400 * 4,
            ],
            "Food": ["apple", "apple", "banana", "banana", "orange", "orange", "pear"],
        }
    ).drop([3])

    df["Datetime"] = to_datetime(df["Timestamp"].apply(lambda t: str(t)), unit="s")
    dfg = df.groupby(Grouper(freq="1D", key="Datetime"))

    # have to sort on index because of unstable sort on values xref GH9212
    result = dfg["Food"].value_counts().sort_index()
    expected = dfg["Food"].apply(Series.value_counts).sort_index()
    expected.index.names = result.index.names

    tm.assert_series_equal(result, expected)
