from collections import OrderedDict
from typing import Any

import click

def index() -> str: ...
def _click_to_tree(
    ctx: click.Context, node: click.Command, ancestors: list[click.Command] | None = None
) -> OrderedDict[str, Any]: ...
