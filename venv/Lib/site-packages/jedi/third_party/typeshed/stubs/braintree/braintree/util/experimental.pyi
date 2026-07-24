from typing import TypeVar

_T = TypeVar("_T")

def Experimental(cls: _T) -> _T: ...
