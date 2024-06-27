from typing import Any, Optional

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


class LoadBalancerNotFoundError(ELBClientError):
    def __init__(self) -> None:
        super().__init__(
            "LoadBalancerNotFound", "The specified load balancer does not exist."
        )


class ListenerNotFoundError(ELBClientError):
    def __init__(self) -> None:
        super().__init__("ListenerNotFound", "The specified listener does not exist.")


class SubnetNotFoundError(ELBClientError):
    def __init__(self) -> None:
        super().__init__("SubnetNotFound", "The specified subnet does not exist.")


class TargetGroupNotFoundError(ELBClientError):
    def __init__(self) -> None:
        super().__init__("TargetGroupNotFound", "One or more target groups not found")


class TooManyTagsError(ELBClientError):
    def __init__(self) -> None:
        super().__init__(
            "TooManyTagsError",
            "The quota for the number of tags that can be assigned to a load balancer has been reached",
        )


class BadHealthCheckDefinition(ELBClientError):
    def __init__(self) -> None:
        super().__init__(
            "ValidationError",
            "HealthCheck Target must begin with one of HTTP, TCP, HTTPS, SSL",
        )


class DuplicateListenerError(ELBClientError):
    def __init__(self) -> None:
        super().__init__(
            "DuplicateListener", "A listener with the specified port already exists."
        )


class DuplicateLoadBalancerName(ELBClientError):
    def __init__(self) -> None:
        super().__init__(
            "DuplicateLoadBalancerName",
            "A load balancer with the specified name already exists.",
        )


class DuplicateTargetGroupName(ELBClientError):
    def __init__(self) -> None:
        super().__init__(
            "DuplicateTargetGroupName",
            "A target group with the specified name already exists.",
        )


class InvalidTargetError(ELBClientError):
    def __init__(self) -> None:
        super().__init__(
            "InvalidTarget",
            "The specified target does not exist or is not in the same VPC as the target group.",
        )


class EmptyListenersError(ELBClientError):
    def __init__(self) -> None:
        super().__init__("ValidationError", "Listeners cannot be empty")


class PriorityInUseError(ELBClientError):
    def __init__(self) -> None:
        super().__init__("PriorityInUse", "The specified priority is in use.")


class InvalidConditionFieldError(ELBClientError):
    VALID_FIELDS = [
        "path-pattern",
        "host-header",
        "http-header",
        "http-request-method",
        "query-string",
        "source-ip",
    ]

    def __init__(self, invalid_name: str):
        valid = ",".join(self.VALID_FIELDS)
        super().__init__(
            "ValidationError",
            f"Condition field '{invalid_name}' must be one of '[{valid}]'",
        )


class InvalidConditionValueError(ELBClientError):
    def __init__(self, msg: str):
        super().__init__("ValidationError", msg)


class InvalidActionTypeError(ELBClientError):
    def __init__(self, invalid_name: Any, index: int):
        super().__init__(
            "ValidationError",
            f"1 validation error detected: Value '{invalid_name}' at 'actions.{index}.member.type' failed to satisfy constraint: Member must satisfy enum value set: [forward, redirect, fixed-response]",
        )


class ActionTargetGroupNotFoundError(ELBClientError):
    def __init__(self, arn: str):
        super().__init__("TargetGroupNotFound", f"Target group '{arn}' not found")


class ListenerOrBalancerMissingError(ELBClientError):
    def __init__(self) -> None:
        super().__init__(
            "ValidationError",
            "You must specify either listener ARNs or a load balancer ARN",
        )


class InvalidDescribeRulesRequest(ELBClientError):
    def __init__(self, msg: str):
        super().__init__("ValidationError", msg)


class ResourceInUseError(ELBClientError):
    def __init__(self, msg: str = "A specified resource is in use"):
        super().__init__("ResourceInUse", msg)


class RuleNotFoundError(ELBClientError):
    def __init__(self, msg: Optional[str] = None):
        msg = msg or "The specified rule does not exist."
        super().__init__("RuleNotFound", msg)


class DuplicatePriorityError(ELBClientError):
    def __init__(self, invalid_value: str):
        super().__init__(
            "ValidationError", f"Priority '{invalid_value}' was provided multiple times"
        )


class InvalidTargetGroupNameError(ELBClientError):
    def __init__(self, msg: str):
        super().__init__("ValidationError", msg)


class InvalidModifyRuleArgumentsError(ELBClientError):
    def __init__(self) -> None:
        super().__init__(
            "ValidationError", "Either conditions or actions must be specified"
        )


class InvalidStatusCodeActionTypeError(ELBClientError):
    def __init__(self, msg: str):
        super().__init__("ValidationError", msg)


class InvalidLoadBalancerActionException(ELBClientError):
    def __init__(self, msg: str):
        super().__init__("InvalidLoadBalancerAction", msg)


class ValidationError(ELBClientError):
    def __init__(self, msg: str):
        super().__init__("ValidationError", msg)


class InvalidConfigurationRequest(ELBClientError):
    def __init__(self, msg: str):
        super().__init__("InvalidConfigurationRequest", msg)
