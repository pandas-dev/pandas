"""
Constants for formatting. Usable by both pandas.core and elsewhere.

The names are chosen to match colorama/ansi.py, whose license is included
in LICENSES/COLORAMA_LICENSE
"""
import re
from typing import Any, Callable, List, Optional

from pandas._libs import missing as libmissing

CSI = "\033["
ANSI_PAT = re.compile(r"\x1B[@-_][0-?]*[ -/]*[@-~]")


def strip_ansi(x):
    return ANSI_PAT.sub("", x)


def format_with(value: str, formatters: List[str]):
    return "".join(formatters + [value, AnsiStyle.RESET_ALL])


class AnsiFore:
    BLACK = f"{CSI}30m"
    RED = f"{CSI}31m"
    GREEN = f"{CSI}32m"
    YELLOW = f"{CSI}33m"
    BLUE = f"{CSI}34m"
    MAGENTA = f"{CSI}35m"
    CYAN = f"{CSI}36m"
    WHITE = f"{CSI}37m"
    RESET = f"{CSI}39m"


class AnsiBack:
    BLACK = f"{CSI}40m"
    RED = f"{CSI}41m"
    GREEN = f"{CSI}42m"
    YELLOW = f"{CSI}43m"
    BLUE = f"{CSI}44m"
    MAGENTA = f"{CSI}45m"
    CYAN = f"{CSI}46m"
    WHITE = f"{CSI}47m"
    RESET = f"{CSI}49m"


class AnsiStyle:
    BRIGHT = f"{CSI}1m"
    DIM = f"{CSI}2m"
    NORMAL = f"{CSI}22m"
    RESET_ALL = f"{CSI}0m"


class NAFormatterMixin:
    def _formatter(
        self, boxed: bool = False, terminal=False
    ) -> Callable[[Any], Optional[str]]:
        def formatter(x):
            if x is libmissing.NA and terminal:
                return format_with("NA", [AnsiFore.RED])
            elif boxed:
                return str(x)
            else:
                return repr(x)

        return formatter
