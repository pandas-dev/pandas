from typing_extensions import TypeAlias

_Graph: TypeAlias = dict[str, list[str | None]]

ADJACENCY_GRAPHS: dict[str, _Graph]
