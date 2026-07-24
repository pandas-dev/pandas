from _typeshed import Incomplete
from typing import Final

from .encryption import StandardSecurityHandler
from .output import ContentWithoutID, OutputProducer
from .syntax import PDFContentStream, PDFObject

HINT_STREAM_OFFSET_LENGTH_PLACEHOLDER: Final[str]
FIRST_PAGE_END_OFFSET_PLACEHOLDER: Final[str]
MAIN_XREF_1ST_ENTRY_OFFSET_PLACEHOLDER: Final[str]
FILE_LENGTH_PLACEHOLDER: Final[str]

class PDFLinearization(PDFObject):
    linearized: str
    n: int
    h: str
    o: Incomplete | None
    e: str
    t: str
    l: str
    def __init__(self, pages_count: int) -> None: ...

class PDFXrefAndTrailer(ContentWithoutID):
    PREV_MAIN_XREF_START_PLACEHOLDER: str
    output_builder: Incomplete
    count: int
    start_obj_id: int
    catalog_obj: Incomplete | None
    info_obj: Incomplete | None
    first_xref: Incomplete | None
    main_xref: Incomplete | None
    startxref: Incomplete | None
    def __init__(self, output_builder) -> None: ...
    @property
    def is_first_xref(self) -> bool: ...
    @property
    def is_main_xref(self) -> bool: ...
    def serialize(self, _security_handler: StandardSecurityHandler | None = None) -> str: ...

class PDFHintStream(PDFContentStream):
    s: Incomplete | None
    t: Incomplete | None
    o: Incomplete | None
    a: Incomplete | None
    e: Incomplete | None
    v: Incomplete | None
    i: Incomplete | None
    c: Incomplete | None
    l: Incomplete | None
    r: Incomplete | None
    b: Incomplete | None

class LinearizedOutputProducer(OutputProducer):
    def bufferize(self) -> bytearray: ...
