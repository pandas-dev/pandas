"""Exceptions raised by the servicediscovery service."""

from moto.core.exceptions import JsonRESTError


class OperationNotFound(JsonRESTError):
    def __init__(self) -> None:
        super().__init__("OperationNotFound", "")


class NamespaceNotFound(JsonRESTError):
    def __init__(self, ns_id: str):
        super().__init__("NamespaceNotFound", f"{ns_id}")


class ServiceNotFound(JsonRESTError):
    def __init__(self, ns_id: str):
        super().__init__("ServiceNotFound", f"{ns_id}")


class ConflictingDomainExists(JsonRESTError):
    def __init__(self, vpc_id: str):
        super().__init__("ConflictingDomainExists", f"{vpc_id}")
