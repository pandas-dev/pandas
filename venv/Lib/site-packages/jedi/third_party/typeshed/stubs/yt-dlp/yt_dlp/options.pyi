import optparse
from collections.abc import Sequence

def parseOpts(
    overrideArguments: Sequence[str] | None = None, ignore_config_files: str = "if_override"
) -> tuple[_YoutubeDLOptionParser, optparse.Values, list[str]]: ...

class _YoutubeDLOptionParser(optparse.OptionParser):
    ALIAS_DEST: str
    ALIAS_TRIGGER_LIMIT: int
    def __init__(self) -> None: ...
    def parse_known_args(
        self, args: Sequence[str] | None = None, values: optparse.Values | None = None, strict: bool = True
    ) -> tuple[optparse.Values, list[str]]: ...

def create_parser() -> _YoutubeDLOptionParser: ...
