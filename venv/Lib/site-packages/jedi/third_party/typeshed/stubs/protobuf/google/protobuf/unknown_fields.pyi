from typing import Any, final

from google.protobuf.message import Message

@final
class UnknownFieldSet:
    def __new__(cls, msg: Message) -> UnknownFieldSet: ...  # noqa: Y034
    def __getitem__(self, index: int, /) -> Any: ...  # Any: internal unknown field object
    def __len__(self) -> int: ...
