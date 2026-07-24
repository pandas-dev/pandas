from moto.core.exceptions import ServiceException


class ECSException(ServiceException):
    pass


class ServiceNotFoundException(ECSException):
    def __init__(self) -> None:
        super().__init__("ServiceNotFoundException", "Service not found.")


class TaskDefinitionNotFoundException(ECSException):
    def __init__(self) -> None:
        super().__init__(
            "ClientException",
            "Unable to describe task definition.",
        )


class RevisionNotFoundException(ECSException):
    def __init__(self) -> None:
        super().__init__("ClientException", "Revision is missing.")


class TaskSetNotFoundException(ECSException):
    def __init__(self) -> None:
        super().__init__(
            "ClientException",
            "The specified task set does not exist.",
        )


class ClusterNotFoundException(ECSException):
    def __init__(self) -> None:
        super().__init__("ClusterNotFoundException", "Cluster not found.")


class EcsClientException(ECSException):
    def __init__(self, message: str):
        super().__init__("ClientException", message)


class InvalidParameterException(ECSException):
    code = "InvalidParameterException"

    def __init__(self, message: str):
        super().__init__(message)


class UnknownAccountSettingException(InvalidParameterException):
    def __init__(self) -> None:
        super().__init__(
            "unknown should be one of [serviceLongArnFormat,taskLongArnFormat,containerInstanceLongArnFormat,containerLongArnFormat,awsvpcTrunking,containerInsights,dualStackIPv6]"
        )


class TaskDefinitionMemoryError(ECSException):
    def __init__(self, container_name: str) -> None:
        super().__init__(
            "ClientException",
            f"Invalid setting for container '{container_name}'. At least one of 'memory' or 'memoryReservation' must be specified.",
        )


class TaskDefinitionMissingPropertyError(ECSException):
    def __init__(self, missing_prop: str) -> None:
        super().__init__(
            "ClientException",
            f"Container.{missing_prop} should not be null or empty.",
        )
