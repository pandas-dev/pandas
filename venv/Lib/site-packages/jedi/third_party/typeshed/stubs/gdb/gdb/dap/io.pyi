from _typeshed import SupportsFlush, SupportsRead, SupportsReadline, SupportsWrite
from typing import Any, type_check_only

import gdb

from .server import _JSONValue
from .startup import DAPQueue

@type_check_only
class _SupportsReadAndReadlineBytes(SupportsRead[bytes], SupportsReadline[bytes]): ...

@type_check_only
class _SupportsWriteAndFlushBytes(SupportsWrite[bytes], SupportsFlush): ...

def read_json(stream: _SupportsReadAndReadlineBytes) -> Any: ...  # returns result of json.loads
def start_json_writer(stream: _SupportsWriteAndFlushBytes, queue: DAPQueue[_JSONValue]) -> gdb.Thread: ...
