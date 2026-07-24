import requests

class SigV4Auth:
    access_key: str
    secret_key: str
    session_token: str | None
    region: str
    def __init__(self, access_key: str, secret_key: str, session_token: str | None = None, region: str = "us-east-1") -> None: ...
    def add_auth(self, request: requests.PreparedRequest) -> None: ...

def generate_sigv4_auth_request(header_value: str | None = None): ...
