from collections.abc import Callable, Iterable, Mapping
from typing import TypeAlias

import optype.numpy as onp

from ._typing import Fetcher

_Datasets: TypeAlias = Fetcher | list[Callable[[], onp.ArrayND]] | tuple[Fetcher, ...]

# NOTE: The implementation explicitly checks that `datasets` sequence is either a `list` or `tuple`, and will raise an
# `AssertionError` for any other (sequence) types.
def _clear_cache(
    datasets: _Datasets | None, cache_dir: str | None = None, method_map: Mapping[str, Iterable[str]] | None = None
) -> None: ...  # undocumented
def clear_cache(datasets: _Datasets | None = None) -> None: ...
