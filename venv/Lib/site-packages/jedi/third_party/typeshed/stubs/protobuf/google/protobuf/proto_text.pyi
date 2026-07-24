from collections.abc import Callable
from typing import TypeVar
from typing_extensions import TypeAlias

from google.protobuf.descriptor_pool import DescriptorPool
from google.protobuf.message import Message

_MessageT = TypeVar("_MessageT", bound=Message)
_MsgFormatter: TypeAlias = Callable[[Message, int, bool], str | None]

def serialize(
    message: Message,
    as_utf8: bool = True,
    as_one_line: bool = False,
    use_short_repeated_primitives: bool = False,
    pointy_brackets: bool = False,
    use_index_order: bool = False,
    use_field_number: bool = False,
    descriptor_pool: DescriptorPool | None = None,
    indent: int = 0,
    message_formatter: _MsgFormatter | None = None,
    print_unknown_fields: bool = False,
    force_colon: bool = False,
) -> str: ...
def parse(
    message_class: type[_MessageT],
    text: str | bytes,
    allow_unknown_extension: bool = False,
    allow_field_number: bool = False,
    descriptor_pool: DescriptorPool | None = None,
    allow_unknown_field: bool = False,
) -> _MessageT: ...
