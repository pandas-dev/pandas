"""Exceptions raised by the servicequotas service."""

from moto.core.exceptions import JsonRESTError


class NoSuchResource(JsonRESTError):
    def __init__(self) -> None:
        super().__init__(
            "NoSuchResourceException",
            "This service is not available in the current Region. Choose a different Region or a different service.",
        )
