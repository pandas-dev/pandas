"""Track current scope to easily calculate the corresponding fine-grained target.

TODO: Use everywhere where we track targets, including in mypy.errors.
"""

from contextlib import contextmanager
from typing import List, Optional, Iterator, Tuple

from mypy.nodes import TypeInfo, FuncBase


SavedScope = Tuple[str, Optional[TypeInfo], Optional[FuncBase]]


class Scope:
    """Track which target we are processing at any given time."""

    def __init__(self) -> None:
        self.module = None  # type: Optional[str]
        self.classes = []  # type: List[TypeInfo]
        self.function = None  # type: Optional[FuncBase]
        # Number of nested scopes ignored (that don't get their own separate targets)
        self.ignored = 0

    def current_module_id(self) -> str:
        assert self.module
        return self.module

    def current_target(self) -> str:
        """Return the current target (non-class; for a class return enclosing module)."""
        assert self.module
        if self.function:
            fullname = self.function.fullname
            return fullname or ''
        return self.module

    def current_full_target(self) -> str:
        """Return the current target (may be a class)."""
        assert self.module
        if self.function:
            return self.function.fullname
        if self.classes:
            return self.classes[-1].fullname
        return self.module

    def current_type_name(self) -> Optional[str]:
        """Return the current type's short name if it exists"""
        return self.classes[-1].name if self.classes else None

    def current_function_name(self) -> Optional[str]:
        """Return the current function's short name if it exists"""
        return self.function.name if self.function else None

    def enter_file(self, prefix: str) -> None:
        self.module = prefix
        self.classes = []
        self.function = None
        self.ignored = 0

    def enter_function(self, fdef: FuncBase) -> None:
        if not self.function:
            self.function = fdef
        else:
            # Nested functions are part of the topmost function target.
            self.ignored += 1

    def enter_class(self, info: TypeInfo) -> None:
        """Enter a class target scope."""
        if not self.function:
            self.classes.append(info)
        else:
            # Classes within functions are part of the enclosing function target.
            self.ignored += 1

    def leave(self) -> None:
        """Leave the innermost scope (can be any kind of scope)."""
        if self.ignored:
            # Leave a scope that's included in the enclosing target.
            self.ignored -= 1
        elif self.function:
            # Function is always the innermost target.
            self.function = None
        elif self.classes:
            # Leave the innermost class.
            self.classes.pop()
        else:
            # Leave module.
            assert self.module
            self.module = None

    def save(self) -> SavedScope:
        """Produce a saved scope that can be entered with saved_scope()"""
        assert self.module
        # We only save the innermost class, which is sufficient since
        # the rest are only needed for when classes are left.
        cls = self.classes[-1] if self.classes else None
        return (self.module, cls, self.function)

    @contextmanager
    def function_scope(self, fdef: FuncBase) -> Iterator[None]:
        self.enter_function(fdef)
        yield
        self.leave()

    @contextmanager
    def class_scope(self, info: TypeInfo) -> Iterator[None]:
        self.enter_class(info)
        yield
        self.leave()

    @contextmanager
    def saved_scope(self, saved: SavedScope) -> Iterator[None]:
        module, info, function = saved
        self.enter_file(module)
        if info:
            self.enter_class(info)
        if function:
            self.enter_function(function)
        yield
        if function:
            self.leave()
        if info:
            self.leave()
        self.leave()
