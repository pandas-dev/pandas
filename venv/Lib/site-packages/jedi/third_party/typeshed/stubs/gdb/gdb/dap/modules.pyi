from _typeshed import Unused
from typing import TypedDict, type_check_only
from typing_extensions import NotRequired

import gdb

@type_check_only
class _Module(TypedDict):
    id: str | None
    name: str | None
    path: NotRequired[str | None]

@type_check_only
class _ModulesResult(TypedDict):
    modules: list[_Module]
    totalModules: int

def module_id(objfile: gdb.Objfile) -> str | None: ...
def is_module(objfile: gdb.Objfile) -> bool: ...
def make_module(objf: gdb.Objfile) -> _Module: ...
def modules(*, startModule: int = 0, moduleCount: int = 0, **args: Unused) -> _ModulesResult: ...
