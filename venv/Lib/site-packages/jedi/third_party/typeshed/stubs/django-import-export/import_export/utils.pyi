from types import TracebackType
from typing import Any, overload

from django.db.models import Field as DjangoField, ForeignObjectRel, Model
from django.db.transaction import Atomic

class atomic_if_using_transaction:
    using_transactions: bool
    context_manager: Atomic
    def __init__(self, using_transactions: bool, using: str | None) -> None: ...
    def __enter__(self) -> None: ...
    def __exit__(
        self, exc_type: type[BaseException] | None, exc_value: BaseException | None, exc_tb: TracebackType | None
    ) -> None: ...

@overload
def get_related_model(field: ForeignObjectRel) -> Model: ...
@overload
def get_related_model(field: DjangoField[Any, Any]) -> Model | None: ...
