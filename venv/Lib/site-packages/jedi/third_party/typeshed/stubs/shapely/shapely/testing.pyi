from ._typing import ArrayLike, OptGeoArrayLike

__all__ = ["assert_geometries_equal"]

def assert_geometries_equal(
    x: OptGeoArrayLike,
    y: OptGeoArrayLike,
    tolerance: ArrayLike[float] = 1e-7,
    equal_none: bool = True,
    equal_nan: bool = True,
    normalize: bool = False,
    err_msg: str = "",
    verbose: bool = True,
) -> None: ...
