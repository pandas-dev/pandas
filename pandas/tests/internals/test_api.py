"""
Tests for the pseudo-public API implemented in internals/api.py and exposed
in core.internals
"""

import pandas as pd
import pandas._testing as tm
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
        "array_manager",
        "base",
        "api",
        "ops",
    ]
    expected = [
        "Block",
        "DatetimeTZBlock",
        "ExtensionBlock",
        "make_block",
        "DataManager",
        "ArrayManager",
        "BlockManager",
        "SingleDataManager",
        "SingleBlockManager",
        "SingleArrayManager",
        "concatenate_managers",
    ]

    result = [x for x in dir(internals) if not x.startswith("__")]
    assert set(result) == set(expected + modules)


def test_make_block_2d_with_dti():
    # GH#41168
    dti = pd.date_range("2012", periods=3, tz="UTC")
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
    with tm.assert_produces_warning(DeprecationWarning, match=msg):
        internals.create_block_manager_from_blocks
