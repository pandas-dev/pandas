import types
from _typeshed import Unused

###

__all__: list[str] = []

def _check_scalar(actual: object, desired: object, *, xp: types.ModuleType | None = None, **kwargs: Unused) -> None: ...
def xp_assert_equal(actual: object, desired: object, *, check_0d: bool = False, xp: types.ModuleType | None = None) -> None: ...
def xp_assert_close(
    actual: object,
    desired: object,
    *,
    check_0d: bool = False,
    rtol: float | None = None,
    atol: float = 0,
    check_namespace: bool = True,
    check_dtype: bool = True,
    check_shape: bool = True,
    err_msg: str = "",
    xp: types.ModuleType | None = None,
) -> None: ...
def xp_assert_less(
    actual: object,
    desired: object,
    *,
    check_0d: bool = False,
    check_namespace: bool = True,
    err_msg: str = "",
    verbose: bool = True,
    xp: types.ModuleType | None = None,
) -> None: ...
def assert_array_almost_equal(
    actual: object,
    desired: object,
    decimal: int = 6,
    *,
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
    *,
    check_namespace: bool = True,
    check_0d: bool = False,
    err_msg: str = "",
    verbose: bool = True,
    xp: types.ModuleType | None = None,
) -> None: ...
