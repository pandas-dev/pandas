from collections.abc import Callable
from typing import Final

# Gleaned from https://github.com/google/jsonnet/blob/master/python/_jsonnet.c
version: Final[str]

def evaluate_file(
    filename: str,
    jpathdir: str | list[str] | None = ...,
    max_stack: int = 500,
    gc_min_objects: int = 1000,
    gc_growth_trigger: float = 2,
    ext_vars: dict[str, str] | None = ...,
    ext_codes: dict[str, str] | None = ...,
    tla_vars: dict[str, str] | None = ...,
    tla_codes: dict[str, str] | None = ...,
    max_trace: int = 20,
    import_callback: Callable[[str, str], tuple[str, object | None]] = ...,
    native_callbacks: dict[str, tuple[tuple[str, ...], Callable[..., object]]] | None = ...,
) -> str: ...
def evaluate_snippet(
    filename: str,
    src: str,
    jpathdir: str | list[str] | None = ...,
    max_stack: int = 500,
    gc_min_objects: int = 1000,
    gc_growth_trigger: float = 2,
    ext_vars: dict[str, str] | None = ...,
    ext_codes: dict[str, str] | None = ...,
    tla_vars: dict[str, str] | None = ...,
    tla_codes: dict[str, str] | None = ...,
    max_trace: int = 20,
    import_callback: Callable[[str, str], tuple[str, object | None]] = ...,
    native_callbacks: dict[str, tuple[tuple[str, ...], Callable[..., object]]] | None = ...,
) -> str: ...
