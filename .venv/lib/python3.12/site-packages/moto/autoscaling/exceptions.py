from moto.core.exceptions import ServiceException


class AutoscalingClientError(ServiceException):
    pass


class ResourceContentionError(AutoscalingClientError):
    code = "ResourceContention"
    message = "You already have a pending update to an Auto Scaling resource (for example, a group, instance, or load balancer)."


class ValidationError(AutoscalingClientError):
    code = "ValidationError"


class InvalidInstanceError(ValidationError):
    def __init__(self, instance_id: str):
        super().__init__(f"Instance [{instance_id}] is invalid.")
