from pandas import DataFrame


class CustomDataFrame(DataFrame):
    """
    Extension of DataFrame as described in [guidelines]

    [guidelines]: https://pandas.pydata.org/pandas-docs/stable/development/extending.html#override-constructor-properties  # noqa
    """

    _metadata = ["_custom_metadata"]

    @property
    def _constructor(self):
        return CustomDataFrame


def test_groupby_with_custom_metadata():
    custom_df = CustomDataFrame(
        [[11, 12, 0], [21, 22, 0], [31, 32, 1]], columns=["a", "b", "g"]
    )
    custom_df._custom_metadata = "Custom metadata"

    for _, group_df in custom_df.groupby("g"):
        assert group_df._custom_metadata == "Custom metadata"
