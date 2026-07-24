from collections.abc import Generator, Sequence
from logging import Logger

LOG: Logger

def expand_paths(
    *, paths: Sequence[str], stdin_display_name: str, filename_patterns: Sequence[str], exclude: Sequence[str]
) -> Generator[str]: ...
