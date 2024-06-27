from .markers import Marker
from .specifiers import SpecifierSet

class InvalidRequirement(ValueError): ...

class Requirement:
    name: str
    url: str | None
    extras: set[str]
    specifier: SpecifierSet
    marker: Marker | None
    def __init__(self, requirement_str: str) -> None: ...
    def __hash__(self) -> int: ...
    def __eq__(self, other: object) -> bool: ...
