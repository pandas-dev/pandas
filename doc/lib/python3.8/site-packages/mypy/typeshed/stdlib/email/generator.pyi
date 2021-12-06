from email.message import Message
from email.policy import Policy
from typing import BinaryIO, Optional, TextIO

class Generator:
    def clone(self, fp: TextIO) -> Generator: ...
    def write(self, s: str) -> None: ...
    def __init__(
        self,
        outfp: TextIO,
        mangle_from_: Optional[bool] = ...,
        maxheaderlen: Optional[int] = ...,
        *,
        policy: Optional[Policy] = ...,
    ) -> None: ...
    def flatten(self, msg: Message, unixfrom: bool = ..., linesep: Optional[str] = ...) -> None: ...

class BytesGenerator:
    def clone(self, fp: BinaryIO) -> BytesGenerator: ...
    def write(self, s: str) -> None: ...
    def __init__(
        self,
        outfp: BinaryIO,
        mangle_from_: Optional[bool] = ...,
        maxheaderlen: Optional[int] = ...,
        *,
        policy: Optional[Policy] = ...,
    ) -> None: ...
    def flatten(self, msg: Message, unixfrom: bool = ..., linesep: Optional[str] = ...) -> None: ...

class DecodedGenerator(Generator):
    def __init__(
        self,
        outfp: TextIO,
        mangle_from_: Optional[bool] = ...,
        maxheaderlen: Optional[int] = ...,
        fmt: Optional[str] = ...,
        *,
        policy: Optional[Policy] = ...,
    ) -> None: ...
