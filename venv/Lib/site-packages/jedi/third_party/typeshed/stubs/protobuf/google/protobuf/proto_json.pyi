from typing import Any, TypeVar

from google.protobuf.descriptor_pool import DescriptorPool
from google.protobuf.message import Message

_MessageT = TypeVar("_MessageT", bound=Message)

def serialize(
    message: Message,
    always_print_fields_with_no_presence: bool = False,
    preserving_proto_field_name: bool = False,
    use_integers_for_enums: bool = False,
    descriptor_pool: DescriptorPool | None = None,
) -> dict[str, Any]: ...
def parse(
    message_class: type[_MessageT],
    js_dict: dict[str, Any],
    ignore_unknown_fields: bool = False,
    descriptor_pool: DescriptorPool | None = None,
    max_recursion_depth: int = 100,
) -> _MessageT: ...
