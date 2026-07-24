import email.message
from typing import Any
from typing import Callable
from typing import Dict
from typing import Generator
from typing import List
from typing import Sequence
from typing import TextIO
from typing import Tuple

def must_decode(value: bytes | str) -> str: ...
def parse(fp: TextIO) -> email.message.Message: ...
def get(msg: email.message.Message, header: str) -> str: ...
def get_all(msg: email.message.Message, header: str) -> List[str]: ...

_header_attr_triple = Tuple[str, str, bool]
_header_attrs = Sequence[_header_attr_triple]
HEADER_ATTRS_1_0: _header_attrs
HEADER_ATTRS_1_1: _header_attrs
HEADER_ATTRS_1_2: _header_attrs
HEADER_ATTRS_2_0 = HEADER_ATTRS_1_2
HEADER_ATTRS_2_1: _header_attrs
HEADER_ATTRS_2_2: _header_attrs
HEADER_ATTRS_2_3: _header_attrs
HEADER_ATTRS_2_4: _header_attrs
HEADER_ATTRS: Dict[str, _header_attrs]
_metadata_version = Tuple[int, int]
METADATA_VERSIONS: Sequence[_metadata_version]
MAX_METADATA_VERSION: _metadata_version
MAX_METADATA_VERSION_STR: str

class UnknownMetadataVersion(UserWarning):
    metadata_version: str
    def __init__(self, metadata_version: str) -> None: ...

class NewMetadataVersion(UserWarning):
    metadata_version: str
    def __init__(self, metadata_version: str) -> None: ...

class Distribution:
    metadata_version: str | None
    name: str | None
    version: str | None
    platforms: Sequence[str]
    supported_platforms: Sequence[str]
    summary: str | None
    description: str | None
    keywords: str | None
    home_page: str | None
    download_url: str | None
    author: str | None
    author_email: str | None
    license: str | None
    classifiers: Sequence[str]
    requires: Sequence[str]
    provides: Sequence[str]
    obsoletes: Sequence[str]
    maintainer: str | None
    maintainer_email: str | None
    requires_python: str | None
    requires_external: Sequence[str]
    requires_dist: Sequence[str]
    provides_dist: Sequence[str]
    obsoletes_dist: Sequence[str]
    project_urls: Sequence[str]
    provides_extras: Sequence[str]
    description_content_type: str | None
    dynamic: Sequence[str]
    license_expression: str | None
    license_file: Sequence[str]
    def extractMetadata(self) -> None: ...
    def read(self) -> bytes: ...
    def parse(self, data: bytes) -> None: ...
    def __iter__(self) -> Generator[str, None, None]: ...
    def iterkeys(self) -> Generator[str, None, None]: ...
