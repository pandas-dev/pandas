from _typeshed import Incomplete
from collections.abc import Awaitable, Callable, Iterable, Mapping
from types import TracebackType
from typing import TypeVar

from .models.dummy_entities import DummySegment, DummySubsegment
from .models.segment import Segment, SegmentContextManager
from .models.subsegment import Subsegment, SubsegmentContextManager
from .recorder import AWSXRayRecorder

_T = TypeVar("_T")

class AsyncSegmentContextManager(SegmentContextManager):
    async def __aenter__(self) -> DummySegment | Segment: ...
    async def __aexit__(
        self, exc_type: type[BaseException] | None, exc_val: BaseException | None, exc_tb: TracebackType | None
    ) -> None: ...

class AsyncSubsegmentContextManager(SubsegmentContextManager):
    async def __call__(
        self, wrapped: Callable[..., Awaitable[_T]], instance, args: Iterable[Incomplete], kwargs: Mapping[str, Incomplete]
    ) -> _T: ...
    async def __aenter__(self) -> DummySubsegment | Subsegment | None: ...
    async def __aexit__(
        self, exc_type: type[BaseException] | None, exc_val: BaseException | None, exc_tb: TracebackType | None
    ) -> None: ...

class AsyncAWSXRayRecorder(AWSXRayRecorder):
    def capture_async(self, name: str | None = None) -> AsyncSubsegmentContextManager: ...
    def in_segment_async(
        self, name: str | None = None, *, traceid: str | None = None, parent_id: str | None = None, sampling: bool | None = None
    ) -> AsyncSegmentContextManager: ...
    def in_subsegment_async(self, name: str | None = None, *, namespace: str = "local") -> AsyncSubsegmentContextManager: ...
    async def record_subsegment_async(
        self,
        wrapped: Callable[..., Awaitable[_T]],
        instance,
        args: Iterable[Incomplete],
        kwargs: Mapping[str, Incomplete],
        name: str,
        namespace: str,
        meta_processor: Callable[..., object] | None,
    ) -> _T: ...
