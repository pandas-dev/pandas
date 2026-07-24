"""Exceptions raised by the appconfig service."""

from moto.core.exceptions import JsonRESTError


class AppNotFoundException(JsonRESTError):
    def __init__(self) -> None:
        super().__init__("ResourceNotFoundException", "Application not found")


class ConfigurationProfileNotFound(JsonRESTError):
    def __init__(self) -> None:
        super().__init__("ResourceNotFoundException", "ConfigurationProfile not found")


class ConfigurationVersionNotFound(JsonRESTError):
    def __init__(self) -> None:
        super().__init__(
            "ResourceNotFoundException", "HostedConfigurationVersion not found"
        )
