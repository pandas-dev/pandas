from collections.abc import Sequence
from re import Pattern
from typing import Protocol, type_check_only

import gdb

def register_xmethod_matcher(
    locus: gdb.Objfile | gdb.Progspace | None, matcher: XMethodMatcher, replace: bool = False
) -> None: ...
@type_check_only
class _XMethod(Protocol):
    name: str
    enabled: bool

class XMethod:
    name: str
    enabled: bool

    def __init__(self, name: str) -> None: ...

class XMethodWorker:
    def get_arg_types(self) -> Sequence[gdb.Type] | gdb.Type | None: ...
    def get_result_type(self, *args: gdb.Value) -> gdb.Type: ...
    def __call__(self, *args: gdb.Value) -> object: ...

class XMethodMatcher:
    name: str
    enabled: bool
    methods: list[_XMethod] | None

    def __init__(self, name: str) -> None: ...
    def match(self, class_type: gdb.Type, method_name: str) -> XMethodWorker | Sequence[XMethodWorker]: ...

@type_check_only
class _SimpleWorkerMethod(Protocol):
    def __call__(self, *args: gdb.Value) -> object: ...

class SimpleXMethodMatcher(XMethodMatcher):
    def __init__(
        self,
        name: str,
        class_matcher: str | Pattern[str],
        method_matcher: str | Pattern[str],
        method_function: _SimpleWorkerMethod,
        *arg_types: Sequence[gdb.Type] | gdb.Type | None,
    ) -> None: ...
