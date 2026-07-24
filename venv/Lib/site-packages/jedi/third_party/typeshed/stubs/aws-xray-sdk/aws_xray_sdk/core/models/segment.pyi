from _typeshed import Incomplete
from types import TracebackType
from typing import Final

from ..recorder import AWSXRayRecorder
from ..utils.atomic_counter import AtomicCounter
from .dummy_entities import DummySegment
from .entity import Entity
from .subsegment import Subsegment

ORIGIN_TRACE_HEADER_ATTR_KEY: Final = "_origin_trace_header"

class SegmentContextManager:
    name: str | None
    segment_kwargs: dict[str, str | bool | None]
    recorder: AWSXRayRecorder
    segment: Segment | None
    def __init__(
        self,
        recorder: AWSXRayRecorder,
        name: str | None = None,
        *,
        traceid: str | None = None,
        parent_id: str | None = None,
        sampling: bool | None = None,
    ) -> None: ...
    def __enter__(self) -> DummySegment | Segment: ...
    def __exit__(
        self, exc_type: type[BaseException] | None, exc_val: BaseException | None, exc_tb: TracebackType | None
    ) -> None: ...

class Segment(Entity):
    trace_id: str
    id: str
    in_progress: bool
    sampled: bool
    user: str | None
    ref_counter: AtomicCounter
    parent_id: str
    service: dict[str, str]
    def __init__(
        self,
        name: str,
        entityid: str | None = None,
        traceid: str | None = None,
        parent_id: str | None = None,
        sampled: bool = True,
    ) -> None: ...
    def add_subsegment(self, subsegment: Subsegment) -> None: ...
    def increment(self) -> None: ...
    def decrement_ref_counter(self) -> None: ...
    def ready_to_send(self) -> bool: ...
    def get_total_subsegments_size(self) -> int: ...
    def decrement_subsegments_size(self) -> int: ...
    def remove_subsegment(self, subsegment: Subsegment) -> None: ...
    def set_user(self, user) -> None: ...
    def set_service(self, service_info: dict[str, str]) -> None: ...
    def set_rule_name(self, rule_name: str) -> None: ...
    def to_dict(self) -> dict[str, Incomplete]: ...
