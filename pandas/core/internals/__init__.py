from pandas.core.internals.api import make_block  # pseudo-public version
from pandas.core.internals.array_manager import (
    ArrayManager,
    SingleArrayManager,
)
from pandas.core.internals.base import (
    DataManager,
    SingleDataManager,
)
from pandas.core.internals.blocks import (  # io.pytables, io.packers
    Block,
    DatetimeBlock,
    DatetimeTZBlock,
    ExtensionBlock,
    FloatBlock,
    NumericBlock,
    ObjectBlock,
    TimeDeltaBlock,
)
from pandas.core.internals.concat import concatenate_managers
from pandas.core.internals.managers import (
    BlockManager,
    SingleBlockManager,
    create_block_manager_from_arrays,
    create_block_manager_from_blocks,
)

__all__ = [
    "Block",
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
    # those two are preserved here for downstream compatibility (GH-33892)
    "create_block_manager_from_arrays",
    "create_block_manager_from_blocks",
]


def __getattr__(name: str):
    import warnings

    if name == "CategoricalBlock":
        warnings.warn(
            "CategoricalBlock is deprecated and will be removed in a future version. "
            "Use ExtensionBlock instead.",
            FutureWarning,
            stacklevel=2,
        )
        return ExtensionBlock

    raise AttributeError(f"module 'pandas.core.internals' has no attribute '{name}'")
