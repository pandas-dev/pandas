from moto.core.exceptions import RESTError


class AutoscalingClientError(RESTError):
    code = 400


class ResourceContentionError(RESTError):
    code = 500

    def __init__(self) -> None:
        super().__init__(
            "ResourceContentionError",
            "You already have a pending update to an Auto Scaling resource (for example, a group, instance, or load balancer).",
        )


class InvalidInstanceError(AutoscalingClientError):
    def __init__(self, instance_id: str):
        super().__init__("ValidationError", f"Instance [{instance_id}] is invalid.")


class ValidationError(AutoscalingClientError):
    def __init__(self, message: str):
        super().__init__("ValidationError", message)
