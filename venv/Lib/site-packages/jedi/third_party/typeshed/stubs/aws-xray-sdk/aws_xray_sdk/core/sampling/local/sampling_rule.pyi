from typing import ClassVar, TypedDict, type_check_only
from typing_extensions import NotRequired

from .reservoir import Reservoir

@type_check_only
class _Rule(TypedDict):
    description: NotRequired[str]
    host: NotRequired[str]
    service_name: NotRequired[str]
    http_method: NotRequired[str]
    url_path: NotRequired[str]
    fixed_target: NotRequired[int]
    rate: NotRequired[float]

class SamplingRule:
    FIXED_TARGET: ClassVar[str]
    RATE: ClassVar[str]
    HOST: ClassVar[str]
    METHOD: ClassVar[str]
    PATH: ClassVar[str]
    SERVICE_NAME: ClassVar[str]
    def __init__(self, rule_dict: _Rule, version: int = 2, default: bool = False) -> None: ...
    def applies(self, host: str | None, method: str | None, path: str | None) -> bool: ...
    @property
    def fixed_target(self) -> int | None: ...
    @property
    def rate(self) -> float | None: ...
    @property
    def host(self) -> str | None: ...
    @property
    def method(self) -> str | None: ...
    @property
    def path(self) -> str | None: ...
    @property
    def reservoir(self) -> Reservoir: ...
    @property
    def version(self): ...
