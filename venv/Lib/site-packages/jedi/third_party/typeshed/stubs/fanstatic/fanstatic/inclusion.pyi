from collections.abc import Iterable

from fanstatic.core import Bundle, NeededResources, Resource

def bundle_resources(resources: Iterable[Resource]) -> list[Resource | Bundle]: ...
def rollup_resources(resources: Iterable[Resource]) -> set[Resource]: ...
def sort_resources(resources: Iterable[Resource]) -> list[Resource]: ...

class Inclusion:
    needed: NeededResources
    resources: list[Resource | Bundle]
    def __init__(
        self,
        needed: NeededResources,
        resources: Iterable[Resource] | None = None,
        compile: bool = False,
        bundle: bool = False,
        mode: str | None = None,
        rollup: bool = False,
    ) -> None: ...
    def __len__(self) -> int: ...
    def render(self) -> str: ...
