from collections.abc import Iterable

from google.protobuf.descriptor import Descriptor
from google.protobuf.descriptor_pb2 import FileDescriptorProto
from google.protobuf.descriptor_pool import DescriptorPool
from google.protobuf.message import Message

def GetMessageClass(descriptor: Descriptor) -> type[Message]: ...
def GetMessageClassesForFiles(files: Iterable[str], pool: DescriptorPool) -> dict[str, type[Message]]: ...

class MessageFactory:
    pool: DescriptorPool
    def __init__(self, pool: DescriptorPool | None = None) -> None: ...

def GetMessages(file_protos: Iterable[FileDescriptorProto], pool: DescriptorPool | None = None) -> dict[str, type[Message]]: ...
