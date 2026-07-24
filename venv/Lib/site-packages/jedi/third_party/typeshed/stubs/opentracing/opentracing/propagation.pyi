from typing import Final

class UnsupportedFormatException(Exception): ...
class InvalidCarrierException(Exception): ...
class SpanContextCorruptedException(Exception): ...

class Format:
    BINARY: Final = "binary"
    TEXT_MAP: Final = "text_map"
    HTTP_HEADERS: Final = "http_headers"
