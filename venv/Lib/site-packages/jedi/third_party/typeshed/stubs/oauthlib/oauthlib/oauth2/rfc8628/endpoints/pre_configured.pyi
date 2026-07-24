from collections.abc import Callable

from oauthlib.openid.connect.core.request_validator import RequestValidator

from .device_authorization import DeviceAuthorizationEndpoint

class DeviceApplicationServer(DeviceAuthorizationEndpoint):
    def __init__(
        self,
        request_validator: RequestValidator,
        verification_uri: str,
        interval: int = 5,
        verification_uri_complete: str | None = None,
        user_code_generator: Callable[[None], str] | None = None,
        **kwargs,
    ) -> None: ...
