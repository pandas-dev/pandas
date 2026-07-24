from collections.abc import Callable
from io import BytesIO
from typing import Any

python_version: tuple[int, int]

BYTES_LITERAL: Callable[[str], bytes]
UNICODE_LITERAL: Callable[[str], str]
BYTES_ORD: Callable[[bytes], int]
BYTES_IO: type[BytesIO]

def fprintf(f: Any, fmt: str, *vargs: Any) -> None: ...

EXCEL_TEXT_TYPES: tuple[type[str], type[bytes], type[bytearray]]
REPR = ascii
xrange = range
unicode: Callable[[bytes, str], str]
ensure_unicode: Callable[[str | bytes], str]
unichr = chr
