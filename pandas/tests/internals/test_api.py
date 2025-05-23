"""
Tests for the pseudo-public API implemented in internals/api.py and exposed
in core.internals
"""

import datetime

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm
from pandas.api.internals import create_dataframe_from_blocks
from pandas.core import internals
from pandas.core.internals import api


def test_internals_api():
    assert internals.make_block is api.make_block


def test_namespace():
    # SUBJECT TO CHANGE

    modules = [
        "blocks",
        "concat",
        "managers",
        "construction",
        "api",
        "ops",
    ]
    expected = [
        "make_block",
        "BlockManager",
        "SingleBlockManager",
        "concatenate_managers",
    ]

    result = [x for x in dir(internals) if not x.startswith("__")]
    assert set(result) == set(expected + modules)


@pytest.mark.parametrize(
    "name",
    [
        "Block",
        "ExtensionBlock",
    ],
)
def test_deprecations(name):
    # GH#55139
    msg = f"{name} is deprecated.* Use public APIs instead"
    with tm.assert_produces_warning(FutureWarning, match=msg):
        getattr(internals, name)


def test_make_block_2d_with_dti():
    # GH#41168
    dti = pd.date_range("2012", periods=3, tz="UTC")

    msg = "make_block is deprecated"
    with tm.assert_produces_warning(DeprecationWarning, match=msg):
        blk = api.make_block(dti, placement=[0])

    assert blk.shape == (1, 3)
    assert blk.values.shape == (1, 3)


def test_create_block_manager_from_blocks_deprecated():
    # GH#33892
    # If they must, downstream packages should get this from internals.api,
    #  not internals.
    msg = (
        "create_block_manager_from_blocks is deprecated and will be "
        "removed in a future version. Use public APIs instead"
    )
    with tm.assert_produces_warning(FutureWarning, match=msg):
        internals.create_block_manager_from_blocks


def test_create_dataframe_from_blocks(float_frame):
    block = float_frame._mgr.blocks[0]
    index = float_frame.index.copy()
    columns = float_frame.columns.copy()

    result = create_dataframe_from_blocks(
        [(block.values, block.mgr_locs.as_array)], index=index, columns=columns
    )
    tm.assert_frame_equal(result, float_frame)


def test_create_dataframe_from_blocks_types():
    df = pd.DataFrame(
        {
            "int": list(range(1, 4)),
            "uint": np.arange(3, 6).astype("uint8"),
            "float": [2.0, np.nan, 3.0],
            "bool": np.array([True, False, True]),
            "boolean": pd.array([True, False, None], dtype="boolean"),
            "string": list("abc"),
            "datetime": pd.date_range("20130101", periods=3),
            "datetimetz": pd.date_range("20130101", periods=3).tz_localize(
                "Europe/Brussels"
            ),
            "timedelta": pd.timedelta_range("1 day", periods=3),
            "period": pd.period_range("2012-01-01", periods=3, freq="D"),
            "categorical": pd.Categorical(["a", "b", "a"]),
            "interval": pd.IntervalIndex.from_tuples([(0, 1), (1, 2), (3, 4)]),
        }
    )

    result = create_dataframe_from_blocks(
        [(block.values, block.mgr_locs.as_array) for block in df._mgr.blocks],
        index=df.index,
        columns=df.columns,
    )
    tm.assert_frame_equal(result, df)


def test_create_dataframe_from_blocks_datetimelike():
    # extension dtypes that have an exact matching numpy dtype can also be
    # be passed as a numpy array
    index, columns = pd.RangeIndex(3), pd.Index(["a", "b", "c", "d"])

    block_array1 = np.arange(
        datetime.datetime(2020, 1, 1),
        datetime.datetime(2020, 1, 7),
        step=datetime.timedelta(1),
    ).reshape((2, 3))
    block_array2 = np.arange(
        datetime.timedelta(1), datetime.timedelta(7), step=datetime.timedelta(1)
    ).reshape((2, 3))
    result = create_dataframe_from_blocks(
        [(block_array1, np.array([0, 2])), (block_array2, np.array([1, 3]))],
        index=index,
        columns=columns,
    )
    expected = pd.DataFrame(
        {
            "a": pd.date_range("2020-01-01", periods=3, unit="us"),
            "b": pd.timedelta_range("1 days", periods=3, unit="us"),
            "c": pd.date_range("2020-01-04", periods=3, unit="us"),
            "d": pd.timedelta_range("4 days", periods=3, unit="us"),
        }
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "array",
    [
        pd.date_range("2020-01-01", periods=3),
        pd.date_range("2020-01-01", periods=3, tz="UTC"),
        pd.period_range("2012-01-01", periods=3, freq="D"),
        pd.timedelta_range("1 day", periods=3),
    ],
)
def test_create_dataframe_from_blocks_1dEA(array):
    # ExtensionArrays can be passed as 1D even if stored under the hood as 2D
    df = pd.DataFrame({"a": array})

    block = df._mgr.blocks[0]
    result = create_dataframe_from_blocks(
        [(block.values[0], block.mgr_locs.as_array)], index=df.index, columns=df.columns
    )
    tm.assert_frame_equal(result, df)
