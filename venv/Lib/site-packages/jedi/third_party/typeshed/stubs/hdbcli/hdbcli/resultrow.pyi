from collections.abc import Iterator, Sequence
from typing import Any, overload
from typing_extensions import disjoint_base

@disjoint_base
class ResultRow:
    def __init__(self, *args: Any, **kwargs: Any) -> None: ...
    column_names: tuple[str, ...]
    column_values: tuple[Any, ...]

    def __len__(self) -> int: ...
    @overload
    def __getitem__(self, index: int, /) -> Any: ...
    @overload
    def __getitem__(self, index: slice, /) -> Sequence[Any]: ...
    def __iter__(self) -> Iterator[Any]: ...
    # __next__, __delitem__, __setitem__ are technically defined but lead always to an error
