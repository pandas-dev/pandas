"""Exceptions raised by the resiliencehub service."""

from moto.core.exceptions import JsonRESTError


class ResourceNotFound(JsonRESTError):
    def __init__(self, msg: str):
        super().__init__("ResourceNotFoundException", msg)


class AppNotFound(ResourceNotFound):
    def __init__(self, arn: str):
        super().__init__(f"App not found for appArn {arn}")


class AppVersionNotFound(ResourceNotFound):
    def __init__(self) -> None:
        super().__init__("App Version not found")


class ResiliencyPolicyNotFound(ResourceNotFound):
    def __init__(self, arn: str):
        super().__init__(f"ResiliencyPolicy {arn} not found")


class ValidationException(JsonRESTError):
    def __init__(self, msg: str):
        super().__init__("ValidationException", msg)
