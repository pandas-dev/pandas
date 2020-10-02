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
