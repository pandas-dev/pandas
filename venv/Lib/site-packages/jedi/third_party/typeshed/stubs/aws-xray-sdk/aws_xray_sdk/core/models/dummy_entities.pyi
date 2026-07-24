from .segment import Segment
from .subsegment import Subsegment

class DummySegment(Segment):
    def __init__(self, name: str = "dummy") -> None: ...

class DummySubsegment(Subsegment):
    def __init__(self, segment, name: str = "dummy") -> None: ...
