from _typeshed import Incomplete
from logging import Logger
from typing import Final, Literal
from typing_extensions import Self, TypeAlias

_SampledTrue: TypeAlias = Literal[True, "1", 1]
_SampledFalse: TypeAlias = Literal[False, "0", 0]
_SampledUnknown: TypeAlias = Literal["?"]
_Sampled: TypeAlias = _SampledTrue | _SampledFalse | _SampledUnknown

log: Logger
ROOT: Final = "Root"
PARENT: Final = "Parent"
SAMPLE: Final = "Sampled"
SELF: Final = "Self"
HEADER_DELIMITER: Final = ";"

class TraceHeader:
    def __init__(
        self,
        root: str | None = None,
        parent: str | None = None,
        sampled: _Sampled | None = None,
        data: dict[str, Incomplete] | None = None,
    ) -> None: ...
    @classmethod
    def from_header_str(cls, header: str | None) -> Self: ...
    def to_header_str(self) -> str: ...
    @property
    def root(self) -> str | None: ...
    @property
    def parent(self) -> str | None: ...
    @property
    def sampled(self) -> Literal[1, 0, "?"] | None: ...
    @property
    def data(self) -> dict[str, Incomplete]: ...
