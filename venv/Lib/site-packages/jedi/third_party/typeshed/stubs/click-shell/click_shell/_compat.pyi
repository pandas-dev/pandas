import types
from collections.abc import Callable
from typing import Any, Final

import click

PY2: Final = False

def get_method_type(func: Callable[..., Any], obj: object) -> types.MethodType: ...
def get_choices(cli: click.Command, prog_name: str, args: list[str], incomplete: str) -> list[str]: ...
