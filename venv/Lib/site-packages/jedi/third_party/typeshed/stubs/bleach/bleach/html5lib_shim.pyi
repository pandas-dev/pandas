import re
from codecs import CodecInfo
from collections.abc import Generator, Iterable, Iterator
from typing import Any, Final, Protocol, type_check_only

# We don't re-export any `html5lib` types / values here, because they are not
# really public and may change at any time. This is just a helper module,
# import things directly from `html5lib` instead!
from html5lib import HTMLParser
from html5lib._inputstream import HTMLBinaryInputStream, HTMLUnicodeInputStream
from html5lib._tokenizer import HTMLTokenizer
from html5lib._trie import Trie
from html5lib.serializer import HTMLSerializer
from html5lib.treewalkers.base import TreeWalker

# Is actually webencodings.Encoding
@type_check_only
class _Encoding(Protocol):
    name: str
    codec_info: CodecInfo
    def __init__(self, name: str, codec_info: CodecInfo) -> None: ...

HTML_TAGS: Final[frozenset[str]]
HTML_TAGS_BLOCK_LEVEL: Final[frozenset[str]]
AMP_SPLIT_RE: Final[re.Pattern[str]]
ENTITIES: Final[dict[str, str]]
ENTITIES_TRIE: Final[Trie]
TAG_TOKEN_TYPES: Final[set[int]]
TAG_TOKEN_TYPE_CHARACTERS: Final[int]
TAG_TOKEN_TYPE_END: Final[int]
TAG_TOKEN_TYPE_PARSEERROR: Final[int]
TAG_TOKEN_TYPE_START: Final[int]

class InputStreamWithMemory:
    position = HTMLUnicodeInputStream.position
    reset = HTMLUnicodeInputStream.reset
    def __init__(self, inner_stream: HTMLUnicodeInputStream) -> None: ...
    @property
    def errors(self) -> list[str]: ...
    @property
    def charEncoding(self) -> tuple[_Encoding, str]: ...
    # If inner_stream wasn't a HTMLBinaryInputStream, this will error at runtime
    # Is a property returning a method, simplified:
    changeEncoding = HTMLBinaryInputStream.changeEncoding
    def char(self) -> str: ...
    def charsUntil(self, characters: Iterable[str], opposite: bool = False) -> str: ...
    def unget(self, char: str | None) -> None: ...
    def get_tag(self) -> str: ...
    def start_tag(self) -> None: ...

class BleachHTMLTokenizer(HTMLTokenizer):
    consume_entities: bool
    stream: InputStreamWithMemory  # type: ignore[assignment]
    emitted_last_token: dict[str, Any] | None
    def __init__(self, consume_entities: bool = False, **kwargs: Any) -> None: ...

class BleachHTMLParser(HTMLParser):
    tags: list[str] | None
    strip: bool
    consume_entities: bool
    def __init__(self, tags: Iterable[str] | None, strip: bool, consume_entities: bool, **kwargs: Any) -> None: ...

class BleachHTMLSerializer(HTMLSerializer):
    escape_rcdata: bool
    def escape_base_amp(self, stoken: str) -> Generator[str]: ...
    def serialize(self, treewalker: TreeWalker, encoding: str | None = None) -> Generator[str]: ...  # type: ignore[override]

def convert_entity(value: str) -> str | None: ...
def convert_entities(text: str) -> str: ...
def match_entity(stream: str) -> str | None: ...
def next_possible_entity(text: str) -> Iterator[str]: ...
