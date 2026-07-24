from typing import Final

esc: Final = "\x1b["
codes: Final[dict[str, str]]
dark_colors: Final[list[str]]
light_colors: Final[list[str]]

def reset_color() -> str: ...
def colorize(color_key: str, text: str) -> str: ...
def ansiformat(attr: str, text: str) -> str: ...
