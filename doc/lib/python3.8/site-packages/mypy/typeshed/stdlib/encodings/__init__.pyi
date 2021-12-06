from codecs import CodecInfo
from typing import Any, Optional, Union

class CodecRegistryError(LookupError, SystemError): ...

def normalize_encoding(encoding: Union[str, bytes]) -> str: ...
def search_function(encoding: str) -> Optional[CodecInfo]: ...

# Needed for submodules
def __getattr__(name: str) -> Any: ...  # incomplete
