from collections import OrderedDict
from typing import Any

from .decoder import TomlDecoder
from .encoder import TomlEncoder

class TomlOrderedDecoder(TomlDecoder[OrderedDict[str, Any]]):
    def __init__(self) -> None: ...

class TomlOrderedEncoder(TomlEncoder[OrderedDict[str, Any]]):
    def __init__(self) -> None: ...
