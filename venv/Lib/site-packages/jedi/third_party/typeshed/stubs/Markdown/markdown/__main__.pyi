import optparse
from logging import Logger
from typing import Any

logger: Logger

def parse_options(
    args: list[str] | None = None, values: optparse.Values | None = None
) -> tuple[dict[str, Any], Any]: ...  # first item is opts dict, second item is Values.verbose field
def run() -> None: ...
