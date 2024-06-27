from _typeshed import Incomplete
from collections.abc import Callable
from typing import TypeVar

Fn = TypeVar("Fn", bound=Callable[..., Incomplete])  # noqa: Y001 # Exists at runtime
__all__ = ("parse_configuration", "read_configuration")

def read_configuration(filepath, find_others: bool = False, ignore_option_errors: bool = False): ...
def parse_configuration(distribution, command_options, ignore_option_errors: bool = False): ...
