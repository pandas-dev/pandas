import types
from typing import Never, TypeAlias

_Ignored: TypeAlias = object

###

__all__: list[str] = []

def _check_scalar(actual: object, desired: object, *, xp: types.ModuleType | None = None, **kwargs: _Ignored) -> None: ...
def xp_assert_equal(actual: object, desired: object, *, check_0d: bool = False, xp: types.ModuleType | None = None) -> None: ...
def xp_assert_close(
    actual: object,
    desired: object,
    *,
    rtol: float | None = None,
    atol: float = 0,
    check_namespace: bool = True,
    check_dtype: bool = True,
    check_shape: bool = True,
    check_0d: bool = False,
    err_msg: str = "",
    xp: types.ModuleType | None = None,
) -> None: ...
def xp_assert_less(
    actual: object,
    desired: object,
    *,
    check_namespace: bool = True,
    check_0d: bool = False,
    err_msg: str = "",
    verbose: bool = True,
    xp: types.ModuleType | None = None,
) -> None: ...
def assert_array_almost_equal(
    actual: object,
    desired: object,
    decimal: int = 6,
    *args: Never,  # see https://github.com/scipy/scipy/pull/24134
    check_namespace: bool = True,
    check_0d: bool = False,
    err_msg: str = "",
    verbose: bool = True,
    xp: types.ModuleType | None = None,
) -> None: ...
def assert_almost_equal(
    actual: object,
    desired: object,
    decimal: int = 7,
    *args: Never,  # see https://github.com/scipy/scipy/pull/24134
    check_namespace: bool = True,
    check_0d: bool = False,
    err_msg: str = "",
    verbose: bool = True,
    xp: types.ModuleType | None = None,
) -> None: ...
