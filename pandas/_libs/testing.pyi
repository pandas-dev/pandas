from collections.abc import Mapping

def assert_dict_equal(a: Mapping, b: Mapping, compare_keys: bool = ...) -> bool: ...
def assert_almost_equal(
    a,
    b,
    rtol: float = ...,
    atol: float = ...,
    check_dtype: bool = ...,
    obj=...,
    lobj=...,
    robj=...,
    index_values=...,
) -> bool: ...
