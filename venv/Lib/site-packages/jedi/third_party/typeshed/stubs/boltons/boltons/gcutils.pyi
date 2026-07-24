from typing import TypeVar

_T = TypeVar("_T")

def get_all(type_obj: type[_T], include_subtypes: bool = True) -> list[_T]: ...

class GCToggler:
    postcollect: bool
    def __init__(self, postcollect: bool = False) -> None: ...
    def __enter__(self) -> None: ...
    def __exit__(self, exc_type, exc_val, exc_tb) -> None: ...

toggle_gc: GCToggler
toggle_gc_postcollect: GCToggler

__all__ = ["get_all", "GCToggler", "toggle_gc", "toggle_gc_postcollect"]
