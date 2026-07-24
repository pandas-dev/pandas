from abc import ABC, abstractmethod

# This interface is unused. Nothing implements `id` as property either.
# IDs can be ints, strings, None, or a Descriptor returning those throughout the codebase.
class ISerialisableFile(ABC):
    @property
    @abstractmethod
    def id(self) -> str | int | None: ...
