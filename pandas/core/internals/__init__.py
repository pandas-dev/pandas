from pandas.core.internals.api import make_block  # 2023-09-18 pyarrow uses this
from pandas.core.internals.concat import concatenate_managers
from pandas.core.internals.managers import (
    BlockManager,
    SingleBlockManager,
)

__all__ = [
    "Block",
    "DatetimeTZBlock",
    "ExtensionBlock",
    "make_block",
    "BlockManager",
    "SingleBlockManager",
    "concatenate_managers",
]


def __getattr__(name: str):
    # GH#55139

    raise AttributeError(f"module 'pandas.core.internals' has no attribute '{name}'")
