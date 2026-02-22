from collections.abc import Callable, Iterable, Sequence
from typing import ClassVar, Final, Literal

import numpy as np
import optype.numpy as onp

__all__ = ["IS_MUSL", "PytestTester", "_TestPythranFunc", "check_free_memory"]

IS_MUSL: Final[bool]
IS_EDITABLE: Final[bool]

class FPUModeChangeWarning(RuntimeWarning): ...

class PytestTester:
    module_name: Final[str]
    def __init__(self, /, module_name: str) -> None: ...
    def __call__(
        self,
        /,
        label: Literal["fast", "full"] = "fast",
        verbose: int = 1,
        extra_argv: Iterable[str] | None = None,
        doctests: bool = False,
        coverage: bool = False,
        tests: Iterable[str] | None = None,
        parallel: int | None = None,
    ) -> bool: ...

class _TestPythranFunc:
    ALL_INTEGER: ClassVar[list[np.int8 | np.int16 | np.int32 | np.int64 | np.intc | np.intp]]
    ALL_FLOAT: ClassVar[list[np.float32 | np.float64]]
    ALL_COMPLEX: ClassVar[list[np.complex64 | np.complex128]]
    arguments: dict[int, tuple[onp.Array, list[np.generic]]]
    partialfunc: Callable[..., object] | None
    expected: object | None
    def setup_method(self, /) -> None: ...
    def get_optional_args(self, /, func: Callable[..., object]) -> dict[str, object]: ...
    def get_max_dtype_list_length(self, /) -> int: ...
    def get_dtype(self, /, dtype_list: Sequence[np.generic], dtype_idx: int) -> np.generic: ...
    def test_all_dtypes(self, /) -> None: ...
    def test_views(self, /) -> None: ...
    def test_strided(self, /) -> None: ...

def check_free_memory(free_mb: float) -> None: ...
