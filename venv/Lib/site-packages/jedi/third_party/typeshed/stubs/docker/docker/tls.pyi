from docker import APIClient

class TLSConfig:
    cert: tuple[str, str] | None
    ca_cert: str | None
    verify: bool | str | None
    def __init__(
        self, client_cert: tuple[str, str] | None = None, ca_cert: str | None = None, verify: bool | str | None = None
    ) -> None: ...
    def configure_client(self, client: APIClient) -> None: ...
