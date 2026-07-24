from typing import Final

from .segment import Segment

MUTATION_UNSUPPORTED_MESSAGE: Final = "FacadeSegments cannot be mutated."

class FacadeSegment(Segment):
    initializing: bool
    def __init__(self, name: str, entityid: str | None, traceid: str | None, sampled: bool | None) -> None: ...
