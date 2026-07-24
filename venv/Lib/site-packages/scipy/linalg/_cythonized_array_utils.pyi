from numpy.typing import NDArray
from typing import Any

def issymmetric(
    a: NDArray[Any],
    atol: None | float = ...,
    rtol: None | float = ...,
) -> bool: ...

def ishermitian(
    a: NDArray[Any],
    atol: None | float = ...,
    rtol: None | float = ...,
) -> bool: ...
