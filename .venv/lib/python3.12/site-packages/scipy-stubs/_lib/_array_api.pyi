import dataclasses
import re
from _typeshed import Incomplete
from collections.abc import Callable, Sequence
from contextlib import _GeneratorContextManager
from types import ModuleType
from typing import Any, Final, Literal, Never, TypeAlias

import _pytest.mark
from array_api_compat import (
    is_array_api_strict_namespace as is_array_api_strict,
    is_cupy_namespace as is_cupy,
    is_jax_namespace as is_jax,
    is_numpy_namespace as is_numpy,
    is_torch_namespace as is_torch,
    numpy as np_compat,
    size as xp_size,
)

_CapabilitiesTable: TypeAlias = dict[Callable[..., Any], dict[str, Any]]

###

__all__ = [
    "SCIPY_ARRAY_API",
    "SCIPY_DEVICE",
    "_asarray",
    "array_namespace",
    "assert_almost_equal",
    "assert_array_almost_equal",
    "default_xp",
    "eager_warns",
    "get_native_namespace_name",
    "is_array_api_strict",
    "is_complex",
    "is_cupy",
    "is_jax",
    "is_lazy_array",
    "is_marray",
    "is_numpy",
    "is_torch",
    "make_xp_pytest_marks",
    "make_xp_pytest_param",
    "make_xp_test_case",
    "np_compat",
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
def default_xp(xp: ModuleType) -> _GeneratorContextManager[None]: ...
def eager_warns(
    warning_type: type[Warning] | tuple[type[Warning], ...], *, match: str | re.Pattern[str] | None = None, xp: ModuleType
) -> _GeneratorContextManager[None]: ...  # _pytest.recwarn.WarningsChecker
def xp_assert_equal(
    actual: object,
    desired: object,
    *,
    check_namespace: bool = True,
    check_dtype: bool = True,
    check_shape: bool = True,
    check_0d: bool = True,
    err_msg: str = "",
    xp: ModuleType | None = None,
) -> None: ...
def xp_assert_close(
    actual: object,
    desired: object,
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
    actual: object,
    desired: object,
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
    actual: object,
    desired: object,
    decimal: int = 6,
    *args: Never,
    check_namespace: bool = True,
    check_0d: bool = True,
    err_msg: str = "",
    xp: ModuleType | None = None,
) -> None: ...
def assert_almost_equal(
    actual: object,
    desired: object,
    decimal: int = 7,
    *args: Never,
    check_namespace: bool = True,
    check_0d: bool = True,
    err_msg: str = "",
    xp: ModuleType | None = None,
) -> None: ...
def xp_unsupported_param_msg(param: Incomplete) -> str: ...
def is_complex(x: Array, xp: ModuleType) -> bool: ...
def get_native_namespace_name(xp: ModuleType) -> str: ...
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
    capabilities_table: _CapabilitiesTable | None = None,
    skip_backends: Sequence[tuple[str, str]] = (),
    xfail_backends: Sequence[tuple[str, str]] = (),
    cpu_only: bool = False,
    np_only: bool = False,
    reason: str | None = None,
    out_of_scope: bool = False,
    exceptions: Sequence[str] = (),
    warnings: Sequence[tuple[str, str]] = (),
    allow_dask_compute: bool = False,
    jax_jit: bool = True,
    extra_note: str | None = None,
) -> dict[str, _XPSphinxCapability]: ...

#
def make_xp_test_case(
    *funcs: Callable[..., Any], capabilities_table: _CapabilitiesTable | None = None
) -> Callable[[Callable[..., None]], Callable[..., None]]: ...
def make_xp_pytest_param(
    func: Callable[..., Any], *args: Any, capabilities_table: _CapabilitiesTable | None = None
) -> _pytest.mark.ParameterSet: ...
def make_xp_pytest_marks(
    *funcs: Callable[..., Any], capabilities_table: _CapabilitiesTable | None = None
) -> Sequence[_pytest.mark.Mark]: ...

xp_capabilities_table: Final[_CapabilitiesTable] = ...

def xp_device_type(a: Array) -> Literal["cpu", "cuda"] | None: ...

# we can't import these from `array_api_compat` due to incomplete annotations, which pyright will complain about
def xp_device(x: object, /) -> Incomplete: ...
def is_lazy_array(x: object) -> bool: ...
