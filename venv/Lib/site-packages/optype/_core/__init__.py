# ruff: noqa: F403
from optype._core import _can, _do, _does, _has, _just
from optype._core._can import *
from optype._core._do import *
from optype._core._does import *
from optype._core._has import *
from optype._core._just import *

__all__: list[str] = []
__all__ += _just.__all__
__all__ += _can.__all__
__all__ += _has.__all__
__all__ += _does.__all__
__all__ += _do.__all__


def __dir__() -> list[str]:
    return __all__


def __rewrite_module_names() -> None:
    name_self = __dir__.__module__
    name_base = name_self.split(".")[0]
    assert name_base != name_self
    assert name_base.isidentifier()

    for module in [_can, _do, _does, _has, _just]:
        for name in module.__all__:
            member: type = getattr(module, name)
            if member.__module__.startswith(name_self):
                member.__module__ = name_base


__rewrite_module_names()
