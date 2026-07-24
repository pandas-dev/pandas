from collections import OrderedDict

from google.protobuf.descriptor_pool import DescriptorPool
from google.protobuf.message import Message

def MakeSimpleProtoClass(
    fields: dict[str, int] | OrderedDict[str, int], full_name: str | None = None, pool: DescriptorPool | None = None
) -> type[Message]: ...
