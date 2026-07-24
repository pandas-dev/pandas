from collections.abc import MutableMapping
from logging import Logger
from typing import Any, Final, overload

log: Logger
SERVICE_NAME: Final = "ec2"
ORIGIN: Final = "AWS::EC2::Instance"
IMDS_URL: Final = "http://169.254.169.254/latest/"

def initialize() -> None: ...
def get_token() -> str | None: ...
def get_metadata(token: str | None = None) -> dict[str, Any]: ...  # result of parse_metadata_json()
def parse_metadata_json(json_str: str | bytes | bytearray) -> dict[str, Any]: ...  # result of json.loads()
@overload
def do_request(url: str, headers: MutableMapping[str, str] | None = None, method: str = "GET") -> str: ...
@overload
def do_request(url: None, headers: MutableMapping[str, str] | None = None, method: str = "GET") -> None: ...
