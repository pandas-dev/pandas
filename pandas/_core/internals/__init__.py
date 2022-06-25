from pandas._core.internals.api import make_block
from pandas._core.internals.array_manager import (
    ArrayManager,
    SingleArrayManager,
)
from pandas._core.internals.base import (
    DataManager,
    SingleDataManager,
)
from pandas._core.internals.blocks import (  # io.pytables, io.packers
    Block,
    DatetimeTZBlock,
    ExtensionBlock,
    NumericBlock,
    ObjectBlock,
)
from pandas._core.internals.concat import concatenate_managers
from pandas._core.internals.managers import (
    BlockManager,
    SingleBlockManager,
    create_block_manager_from_blocks,
)

__all__ = [
    "Block",
    "CategoricalBlock",
    "NumericBlock",
    "DatetimeTZBlock",
    "ExtensionBlock",
    "ObjectBlock",
    "make_block",
    "DataManager",
    "ArrayManager",
    "BlockManager",
    "SingleDataManager",
    "SingleBlockManager",
    "SingleArrayManager",
    "concatenate_managers",
    # this is preserved here for downstream compatibility (GH-33892)
    "create_block_manager_from_blocks",
]


def __getattr__(name: str):
    import warnings

    from pandas.util._exceptions import find_stack_level

    if name == "CategoricalBlock":
        warnings.warn(
            "CategoricalBlock is deprecated and will be removed in a future version. "
            "Use ExtensionBlock instead.",
            DeprecationWarning,
            stacklevel=find_stack_level(),
        )
        from pandas._core.internals.blocks import CategoricalBlock

        return CategoricalBlock

    raise AttributeError(f"module 'pandas._core.internals' has no attribute '{name}'")
