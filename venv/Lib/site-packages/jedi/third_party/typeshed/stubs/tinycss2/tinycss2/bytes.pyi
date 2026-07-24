from webencodings import Encoding

from .ast import Node

def decode_stylesheet_bytes(
    css_bytes: bytes, protocol_encoding: str | None = None, environment_encoding: Encoding | None = None
) -> tuple[str, Encoding]: ...
def parse_stylesheet_bytes(
    css_bytes: bytes,
    protocol_encoding: str | None = None,
    environment_encoding: Encoding | None = None,
    skip_comments: bool = False,
    skip_whitespace: bool = False,
) -> tuple[list[Node], Encoding]: ...
