import logging
from _typeshed import Incomplete
from collections.abc import Callable, Sequence
from contextlib import _GeneratorContextManager
from typing import Any, TypeVar, overload

from ..std import tqdm as std_tqdm

_TqdmT = TypeVar("_TqdmT", bound=std_tqdm[Any])

def logging_redirect_tqdm(
    loggers: Sequence[logging.Logger] | None = None, tqdm_class: type[std_tqdm[Any]] = ...
) -> _GeneratorContextManager[None]: ...

# TODO: type *args, **kwargs here more precisely
@overload
def tqdm_logging_redirect(*args, tqdm_class: Callable[..., _TqdmT], **kwargs) -> _GeneratorContextManager[_TqdmT]: ...
@overload
def tqdm_logging_redirect(*args, **kwargs) -> _GeneratorContextManager[std_tqdm[Incomplete]]: ...
