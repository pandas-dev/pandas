from _typeshed import Incomplete
from datetime import datetime
from typing import TypedDict, type_check_only
from typing_extensions import Self

@type_check_only
class _ParamsDict(TypedDict, total=False):
    ach_mandate_text: str
    ach_mandate_accepted_at: str

class TransactionUsBankAccountRequest:
    parent: Incomplete
    def __init__(self, parent) -> None: ...
    def ach_mandate_text(self, ach_mandate_text: str) -> Self: ...
    def ach_mandate_accepted_at(self, ach_mandate_accepted_at: str | datetime) -> Self: ...
    def done(self): ...
    def to_param_dict(self) -> _ParamsDict: ...
