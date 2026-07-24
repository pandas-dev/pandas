from .base import Client

class BackendApplicationClient(Client):
    grant_type: str
    def prepare_request_body(
        self,
        body: str = "",
        scope: str | set[object] | tuple[object] | list[object] | None = None,
        include_client_id: bool = False,
        *,
        code_verifier: str | None = None,
        client_id: str | None = None,
        client_secret: str | None = None,
        **kwargs,
    ) -> str: ...
