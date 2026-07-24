from typing import TypeVar

_T = TypeVar("_T")

# See note in stubs/auth0-python/@tests/stubtest_allowlist.txt about _async methods
def asyncify(cls: type[_T]) -> type[_T]: ...
