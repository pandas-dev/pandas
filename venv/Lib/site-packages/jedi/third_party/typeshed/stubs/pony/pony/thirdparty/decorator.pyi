import re
from _typeshed import Incomplete
from collections.abc import Callable
from inspect import getfullargspec as getfullargspec
from typing import Final

__version__: Final[str]
__all__ = ["decorator", "FunctionMaker", "contextmanager"]

def get_init(cls: object) -> Callable[..., None]: ...

DEF: re.Pattern[str]

class FunctionMaker:
    shortsignature: Incomplete
    name: Incomplete
    doc: Incomplete
    module: Incomplete
    annotations: Incomplete
    signature: Incomplete
    dict: Incomplete
    defaults: Incomplete
    def __init__(self, func=None, name=None, signature=None, defaults=None, doc=None, module=None, funcdict=None) -> None: ...
    def update(self, func, **kw) -> None: ...
    def make(self, src_templ, evaldict=None, addsource: bool = False, **attrs): ...
    @classmethod
    def create(cls, obj, body, evaldict, defaults=None, doc=None, module=None, addsource: bool = True, **attrs): ...

def decorator(caller, func=None): ...

contextmanager: Incomplete
