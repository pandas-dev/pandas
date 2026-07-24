from _typeshed import Incomplete
from typing_extensions import TypeAlias

TimeoutType: TypeAlias = float | tuple[float, float]
RequestData: TypeAlias = dict[str, Incomplete] | list[Incomplete]
