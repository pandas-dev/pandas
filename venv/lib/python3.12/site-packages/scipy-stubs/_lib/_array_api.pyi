import dataclasses
import re
from _typeshed import Incomplete
from collections.abc import Generator, Sequence
from contextlib import contextmanager
from types import ModuleType
from typing import Any, Final, Literal, TypeAlias

from array_api_compat import (
    device as xp_device,  # pyright: ignore[reportUnknownVariableType]
    is_array_api_strict_namespace as is_array_api_strict,
    is_cupy_namespace as is_cupy,
    is_jax_namespace as is_jax,
    is_lazy_array as is_lazy_array,
    is_numpy_namespace as is_numpy,
    is_torch_namespace as is_torch,
    size as xp_size,
)

__all__ = [
    "SCIPY_ARRAY_API",
    "SCIPY_DEVICE",
    "_asarray",
    "array_namespace",
    "assert_almost_equal",
    "assert_array_almost_equal",
    "default_xp",
    "eager_warns",
    "is_array_api_strict",
    "is_complex",
    "is_cupy",
    "is_jax",
    "is_lazy_array",
    "is_marray",
    "is_numpy",
    "is_torch",
    "scipy_namespace_for",
    "xp_assert_close",
    "xp_assert_equal",
    "xp_assert_less",
    "xp_capabilities",
    "xp_copy",
    "xp_device",
    "xp_promote",
    "xp_ravel",
    "xp_result_type",
    "xp_size",
    "xp_unsupported_param_msg",
    "xp_vector_norm",
]

SCIPY_ARRAY_API: Final[str | Literal[False]] = ...
SCIPY_DEVICE: Final[str] = ...

Array: TypeAlias = Incomplete

def array_namespace(*arrays: Array) -> ModuleType: ...
def _asarray(
    array: Any,
    dtype: Any = None,
    order: Literal["K", "A", "C", "F"] | None = None,
    copy: bool | None = None,
    *,
    xp: ModuleType | None = None,
    check_finite: bool = False,
    subok: bool = False,
) -> Array: ...
def xp_copy(x: Array, *, xp: ModuleType | None = None) -> Array: ...
@contextmanager
def default_xp(xp: ModuleType) -> Generator[None]: ...
def eager_warns(
    x: Array, warning_type: type[Warning] | tuple[type[Warning], ...], match: str | re.Pattern[str] | None = None
) -> Incomplete: ...  # _pytest.recwarn.WarningsChecker
def xp_assert_equal(
    actual: Incomplete,
    desired: Incomplete,
    *,
    check_namespace: bool = True,
    check_dtype: bool = True,
    check_shape: bool = True,
    check_0d: bool = True,
    err_msg: str = "",
    xp: ModuleType | None = None,
) -> None: ...
def xp_assert_close(
    actual: Incomplete,
    desired: Incomplete,
    *,
    rtol: Incomplete | None = None,
    atol: int = 0,
    check_namespace: bool = True,
    check_dtype: bool = True,
    check_shape: bool = True,
    check_0d: bool = True,
    err_msg: str = "",
    xp: ModuleType | None = None,
) -> None: ...
def xp_assert_less(
    actual: Incomplete,
    desired: Incomplete,
    *,
    check_namespace: bool = True,
    check_dtype: bool = True,
    check_shape: bool = True,
    check_0d: bool = True,
    err_msg: str = "",
    verbose: bool = True,
    xp: ModuleType | None = None,
) -> None: ...
def assert_array_almost_equal(
    actual: Incomplete, desired: Incomplete, decimal: int = 6, *args: Incomplete, **kwds: Incomplete
) -> None: ...
def assert_almost_equal(
    actual: Incomplete, desired: Incomplete, decimal: int = 7, *args: Incomplete, **kwds: Incomplete
) -> None: ...
def xp_unsupported_param_msg(param: Incomplete) -> str: ...
def is_complex(x: Array, xp: ModuleType) -> bool: ...
def scipy_namespace_for(xp: ModuleType) -> ModuleType | None: ...
def xp_vector_norm(
    x: Array,
    /,
    *,
    axis: int | tuple[int, ...] | None = None,
    keepdims: bool = False,
    ord: float = 2,
    xp: ModuleType | None = None,
) -> Array: ...
def xp_ravel(x: Array, /, *, xp: ModuleType | None = None) -> Array: ...
def xp_result_type(*args: Incomplete, force_floating: bool = False, xp: ModuleType) -> type: ...
def xp_promote(*args: Incomplete, broadcast: bool = False, force_floating: bool = False, xp: ModuleType | None) -> Array: ...
def is_marray(xp: ModuleType) -> bool: ...

@dataclasses.dataclass(repr=False)
class _XPSphinxCapability:
    cpu: bool | None
    gpu: bool | None
    warnings: list[str] = ...

    def _render(self, /, value: object) -> str: ...

def xp_capabilities(
    *,
    capabilities_table: Incomplete | None = None,
    skip_backends: Sequence[tuple[str, str]] = (),
    xfail_backends: Sequence[tuple[str, str]] = (),
    cpu_only: bool = False,
    np_only: bool = False,
    reason: str | None = None,
    exceptions: Sequence[str] = (),
    warnings: Sequence[tuple[str, str]] = (),
    allow_dask_compute: bool = False,
    jax_jit: bool = True,
) -> dict[str, _XPSphinxCapability]: ...
