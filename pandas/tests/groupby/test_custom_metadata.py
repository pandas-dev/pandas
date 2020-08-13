"""
Test metadata propagation in groupby

The PandasTable class below is implemented according to the [guidelines],
and as such would expect `__finalize__` to always be called so that the
`_pandastable_metadata` is always populated.

[guidelines]: https://pandas.pydata.org/pandas-docs/stable/development/extending.html#override-constructor-properties  # noqa
"""

import pytest
import pandas as pd
from warnings import warn
from typing import List


_TABLE_METADATA_FIELD_NAME = '_pandastable_metadata'


def _combine_metadata(data: List[str]) -> str:
    """
    A mock implementation for testing
    """
    return '+'.join(data)


class PandasTable(pd.DataFrame):
    """
    A pandas dataframe subclass with associated table metadata.
    """

    _metadata = [_TABLE_METADATA_FIELD_NAME]  # Register metadata fieldnames here

    @property
    def _constructor(self):
        return PandasTable

    def __finalize__(self, other, method=None, **kwargs):
        """
        This method be called after constructor to populate metadata

        The "method" argument is subject to change and must be handled robustly.
        """
        src = [other]  # more logic here in actual implementation
        metadata = _combine_metadata(
            [d.get_metadata() for d in src if isinstance(d, PandasTable)])

        if not metadata:
            warn('__finalize__ unable to combine metadata for method "{method}", '
                 'falling back to DataFrame')
            return pd.DataFrame(self)
        object.__setattr__(self, _TABLE_METADATA_FIELD_NAME, metadata)
        return self

    def get_metadata(self):
        metadata = getattr(self, _TABLE_METADATA_FIELD_NAME, None)
        if metadata is None:
            warn('PandasTable object not correctly initialized: no metadata')
        return metadata

    @staticmethod
    def from_table_data(df: pd.DataFrame, metadata) -> 'PandasTable':
        df = PandasTable(df)
        object.__setattr__(df, _TABLE_METADATA_FIELD_NAME, metadata)
        return df


@pytest.fixture
def dft():
    df = pd.DataFrame([[11, 12, 0], [21, 22, 0], [31, 32, 1]], columns={'a', 'b', 'g'})
    return PandasTable.from_table_data(df, 'My metadata')


def test_initial_metadata(dft):
    assert dft.get_metadata() == 'My metadata'


def test_basic_propagation(dft):
    gg = dft.loc[dft.g == 0, :]
    assert gg.get_metadata() == 'My metadata'


def test_groupby(dft):
    gg = [ab for g, ab in dft.groupby('g')]
    assert gg[0].get_metadata() is not None
