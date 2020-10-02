import numpy as np
import pytest

import pandas as pd


class SubclassedDataFrame2(pd.DataFrame):
    """
    Extension of DataFrame as described in [guidelines]

    [guidelines]: https://pandas.pydata.org/pandas-docs/stable/development/extending.html#override-constructor-properties  # noqa
    """

    # normal properties
    _metadata = ["added_property"]

    @property
    def _constructor(self):
        return SubclassedDataFrame2


@pytest.mark.xfail(
    reason="Missing high-performance implementation of metadata support for groupby. "
    "__finalize__ is not called for grouped dataframes"
)
def test_groupby_with_custom_metadata():
    custom_df = SubclassedDataFrame2(
        [[11, 12, 0], [21, 22, 0], [31, 32, 1]], columns=["a", "b", "g"]
    )
    custom_df.added_property = "hello_pandas"
    grouped = custom_df.groupby("g")
    for _, group_df in grouped:
        assert group_df.added_property == "hello_pandas"


def test_groupby_sum_with_custom_metadata():
    my_data_as_dictionary = {
        "mycategorical": [1, 1, 2],
        "myfloat1": [1.0, 2.0, 3.0],
        "myfloat2": [1.0, 2.0, 3.0],
    }
    sdf = SubclassedDataFrame2(my_data_as_dictionary)
    sdf.added_property = "hello pandas"
    grouped = sdf.groupby("mycategorical")[["myfloat1", "myfloat2"]]
    df = grouped.sum()
    assert df.added_property == "hello pandas"


def test_groupby_apply_performance_with_custom_metadata():
    """
    Check that __finalize__ is not called for each group during .apply
    """
    counter = {"value": 0}

    class SubclassedDataFrame3(pd.DataFrame):
        # normal properties
        _metadata = ["added_property"]

        @property
        def _constructor(self):
            return SubclassedDataFrame3

        def __finalize__(self, *args, **kwargs):
            counter["value"] += 1
            return super().__finalize__(*args, **kwargs)

    N = 10 ** 4
    labels = np.random.randint(0, N // 5, size=N)
    labels2 = np.random.randint(0, 3, size=N)
    df = SubclassedDataFrame3(
        {
            "key": labels,
            "key2": labels2,
            "value1": np.random.randn(N),
            "value2": ["foo", "bar", "baz", "qux"] * (N // 4),
        }
    )
    df.index = pd.CategoricalIndex(df.key)
    df.groupby(level="key").apply(lambda x: 1)
    finalize_call_count = counter["value"]
    assert finalize_call_count < 42  # Random number << len(labels2)
