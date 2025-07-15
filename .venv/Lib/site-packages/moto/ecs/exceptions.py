from moto.core.exceptions import JsonRESTError, RESTError


class ServiceNotFoundException(RESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            error_type="ServiceNotFoundException", message="Service not found."
        )


class TaskDefinitionNotFoundException(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            error_type="ClientException",
            message="Unable to describe task definition.",
        )


class RevisionNotFoundException(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(error_type="ClientException", message="Revision is missing.")


class TaskSetNotFoundException(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            error_type="ClientException",
            message="The specified task set does not exist.",
        )


class ClusterNotFoundException(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            error_type="ClusterNotFoundException", message="Cluster not found."
        )


class EcsClientException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__(error_type="ClientException", message=message)


class InvalidParameterException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__(error_type="InvalidParameterException", message=message)


class UnknownAccountSettingException(InvalidParameterException):
    def __init__(self) -> None:
        super().__init__(
            "unknown should be one of [serviceLongArnFormat,taskLongArnFormat,containerInstanceLongArnFormat,containerLongArnFormat,awsvpcTrunking,containerInsights,dualStackIPv6]"
        )


class TaskDefinitionMemoryError(JsonRESTError):
    def __init__(self, container_name: str) -> None:
        super().__init__(
            error_type="ClientException",
            message=f"Invalid setting for container '{container_name}'. At least one of 'memory' or 'memoryReservation' must be specified.",
        )


class TaskDefinitionMissingPropertyError(JsonRESTError):
    def __init__(self, missing_prop: str) -> None:
        super().__init__(
            error_type="ClientException",
            message=f"Container.{missing_prop} should not be null or empty.",
        )
