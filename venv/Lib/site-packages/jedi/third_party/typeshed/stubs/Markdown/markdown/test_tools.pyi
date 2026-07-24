import unittest
from _typeshed import Unused
from typing import Any

__all__ = ["TestCase", "LegacyTestCase", "Kwargs"]

class TestCase(unittest.TestCase):
    default_kwargs: dict[str, Any]  # taken from source code
    def assertMarkdownRenders(
        self,
        source: str,
        expected: str,
        expected_attrs: dict[str, Any] | None = None,  # values passing to self.assertEqual()
        **kwargs,
    ) -> None: ...
    def dedent(self, text: str) -> str: ...

class recursionlimit:
    limit: int
    old_limit: int
    def __init__(self, limit: int) -> None: ...
    def __enter__(self) -> None: ...
    def __exit__(self, type: Unused, value: Unused, tb: Unused) -> None: ...

class Kwargs(dict[str, Any]): ...

class LegacyTestMeta(type):
    def __new__(cls, name: str, bases: tuple[type, ...], dct: dict[str, Any]): ...  # dct is namespace argument for type.__new__()

class LegacyTestCase(unittest.TestCase, metaclass=LegacyTestMeta): ...
