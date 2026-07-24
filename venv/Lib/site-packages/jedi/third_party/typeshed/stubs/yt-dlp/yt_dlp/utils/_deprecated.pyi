import re
from collections.abc import Mapping
from typing import Any
from typing_extensions import deprecated

@deprecated("Use `list` instead.")
def bytes_to_intlist(bs: bytes) -> list[int]: ...
@deprecated("Use `bytes` instead.")
def intlist_to_bytes(xs: list[int]) -> bytes: ...
@deprecated("Use `yt_dlp.utils.jwt_encode` instead.")
def jwt_encode_hs256(payload_data: Any, key: str, headers: Mapping[str, Any] = {}) -> bytes: ...  # Passed to json.dumps().

compiled_regex_type: type[re.Pattern[Any]]
