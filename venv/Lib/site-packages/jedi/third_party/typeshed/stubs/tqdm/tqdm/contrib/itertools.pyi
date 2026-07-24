from _typeshed import Incomplete
from collections.abc import Generator, Iterable

__all__ = ["product"]

def product(*iterables: Iterable[Incomplete], **tqdm_kwargs) -> Generator[Incomplete]: ...
