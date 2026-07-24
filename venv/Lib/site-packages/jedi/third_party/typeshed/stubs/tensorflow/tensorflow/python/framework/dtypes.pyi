import dataclasses
from _typeshed import Incomplete

@dataclasses.dataclass(frozen=True)
class HandleData:
    shape_inference: Incomplete | None = None
    alias_id: int | None = None
