from typing import (
    Sequence,
    overload,
)

from pandas._typing import (
    AnyArrayLike,
    DataFrame,
    Index,
    Series,
)

# These cannot _really_ be just FrameOps/SeriesOps, as those are
#  mixins to DataFrame/Series. We include those here so that the annotations
#  in the mixin are correct.
from pandas.core.ops.methods import (
    FrameOps,
    SeriesOps,
)

# note: this is a lie to make type checkers happy (they special
# case property). cache_readonly uses attribute names similar to
# property (fget) but it does not provide fset and fdel.
cache_readonly = property

class AxisProperty:
    axis: int
    def __init__(self, axis: int = ..., doc: str = ...) -> None: ...
    @overload
    def __get__(
        self, obj: DataFrame | Series | SeriesOps | FrameOps, type
    ) -> Index: ...
    @overload
    def __get__(self, obj: None, type) -> AxisProperty: ...
    def __set__(
        self,
        obj: DataFrame | Series | SeriesOps | FrameOps,
        value: AnyArrayLike | Sequence,
    ) -> None: ...
