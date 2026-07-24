from typing import Any, Final

from jmespath import parser as parser
from jmespath.visitor import Options as Options

__version__: Final[str]

def compile(expression: str) -> parser.ParsedResult: ...
def search(expression: str, data: Any, options: Options | None = None) -> Any: ...
