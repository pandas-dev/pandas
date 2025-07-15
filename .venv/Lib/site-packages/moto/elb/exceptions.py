from typing import Any

from moto.core.exceptions import RESTError


class ELBClientError(RESTError):
    code = 400

    def __init__(self, error_type: str, message: str):
        super().__init__(error_type, message, template="wrapped_single_error")


class DuplicateTagKeysError(ELBClientError):
    def __init__(self, cidr: Any):
        super().__init__(
            "DuplicateTagKeys", f"Tag key was specified more than once: {cidr}"
        )


class CertificateNotFoundException(ELBClientError):
    def __init__(self) -> None:
        super().__init__(
            "CertificateNotFoundException", "Supplied certificate was not found"
        )


class LoadBalancerNotFoundError(ELBClientError):
    def __init__(self, name: str):
        super().__init__(
            "LoadBalancerNotFound",
            f"The specified load balancer does not exist: {name}",
        )


class NoActiveLoadBalancerFoundError(ELBClientError):
    def __init__(self, name: str):
        super().__init__(
            "LoadBalancerNotFound", f"There is no ACTIVE Load Balancer named '{name}'"
        )


class PolicyNotFoundError(ELBClientError):
    def __init__(self) -> None:
        super().__init__(
            "PolicyNotFound", "There is no policy with name . for load balancer ."
        )


class TooManyTagsError(ELBClientError):
    def __init__(self) -> None:
        super().__init__(
            "LoadBalancerNotFound",
            "The quota for the number of tags that can be assigned to a load balancer has been reached",
        )


class BadHealthCheckDefinition(ELBClientError):
    def __init__(self) -> None:
        super().__init__(
            "ValidationError",
            "HealthCheck Target must begin with one of HTTP, TCP, HTTPS, SSL",
        )


class DuplicateListenerError(ELBClientError):
    def __init__(self, name: str, port: str):
        super().__init__(
            "DuplicateListener",
            f"A listener already exists for {name} with LoadBalancerPort {port}, but with a different InstancePort, Protocol, or SSLCertificateId",
        )


class DuplicateLoadBalancerName(ELBClientError):
    def __init__(self, name: str):
        super().__init__(
            "DuplicateLoadBalancerName",
            f"The specified load balancer name already exists for this account: {name}",
        )


class EmptyListenersError(ELBClientError):
    def __init__(self) -> None:
        super().__init__("ValidationError", "Listeners cannot be empty")


class InvalidSecurityGroupError(ELBClientError):
    def __init__(self) -> None:
        super().__init__(
            "ValidationError",
            "One or more of the specified security groups do not exist.",
        )
