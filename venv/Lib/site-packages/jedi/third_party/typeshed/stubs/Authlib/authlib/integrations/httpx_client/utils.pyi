from _typeshed import Incomplete
from typing import Final

HTTPX_CLIENT_KWARGS: Final[list[str]]

def extract_client_kwargs(kwargs) -> dict[str, Incomplete]: ...
def build_request(url, headers, body, initial_request): ...
