"""
Tests for the pseudo-public API implemented in internals/api.py and exposed
in core.internals
"""

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
        "CategoricalBlock",
        "NumericBlock",
        "DatetimeBlock",
        "DatetimeTZBlock",
        "ExtensionBlock",
        "FloatBlock",
        "ObjectBlock",
        "TimeDeltaBlock",
        "make_block",
        "DataManager",
        "ArrayManager",
        "BlockManager",
        "SingleDataManager",
        "SingleBlockManager",
        "SingleArrayManager",
        "concatenate_managers",
        "create_block_manager_from_arrays",
        "create_block_manager_from_blocks",
    ]

    result = [x for x in dir(internals) if not x.startswith("__")]
    assert set(result) == set(expected + modules)
