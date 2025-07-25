from typing import Any

from moto.core.exceptions import JsonRESTError


class LambdaClientError(JsonRESTError):
    def __init__(self, error: str, message: str):
        super().__init__(error, message)


class CrossAccountNotAllowed(LambdaClientError):
    def __init__(self) -> None:
        super().__init__(
            "AccessDeniedException", "Cross-account pass role is not allowed."
        )


class FunctionAlreadyExists(LambdaClientError):
    code = 409

    def __init__(self, function_name: str) -> None:
        message = f"Function already exist: {function_name}"
        super().__init__("ResourceConflictException", message)


class InvalidParameterValueException(LambdaClientError):
    def __init__(self, message: str):
        super().__init__("InvalidParameterValueException", message)


class InvalidRoleFormat(LambdaClientError):
    pattern = r"arn:(aws[a-zA-Z-]*)?:iam::(\d{12}):role/?[a-zA-Z_0-9+=,.@\-_/]+"

    def __init__(self, role: str):
        message = f"1 validation error detected: Value '{role}' at 'role' failed to satisfy constraint: Member must satisfy regular expression pattern: {InvalidRoleFormat.pattern}"
        super().__init__("ValidationException", message)


class PreconditionFailedException(JsonRESTError):
    code = 412

    def __init__(self, message: str):
        super().__init__("PreconditionFailedException", message)


class ConflictException(LambdaClientError):
    code = 409

    def __init__(self, message: str):
        super().__init__("ConflictException", message)


class UnknownAliasException(LambdaClientError):
    code = 404

    def __init__(self, arn: str):
        super().__init__("ResourceNotFoundException", f"Cannot find alias arn: {arn}")


class UnknownFunctionException(LambdaClientError):
    code = 404

    def __init__(self, arn: str):
        super().__init__("ResourceNotFoundException", f"Function not found: {arn}")


class GenericResourcNotFound(LambdaClientError):
    code = 404

    def __init__(self) -> None:
        super().__init__(
            "ResourceNotFoundException", "The resource you requested does not exist."
        )


class UnknownLayerException(LambdaClientError):
    code = 404

    def __init__(self) -> None:
        super().__init__("ResourceNotFoundException", "Cannot find layer")


class UnknownLayerVersionException(LambdaClientError):
    code = 404

    def __init__(self, arns: Any) -> None:
        super().__init__(
            "ResourceNotFoundException",
            f"One or more LayerVersion does not exist {arns}",
        )


class UnknownPolicyException(LambdaClientError):
    code = 404

    def __init__(self) -> None:
        super().__init__(
            "ResourceNotFoundException",
            "No policy is associated with the given resource.",
        )


class UnknownEventConfig(LambdaClientError):
    code = 404

    def __init__(self, arn: str) -> None:
        super().__init__(
            "ResourceNotFoundException",
            f"The function {arn} doesn't have an EventInvokeConfig",
        )


class ValidationException(LambdaClientError):
    def __init__(self, value: str, property_name: str, specific_message: str):
        message = f"1 validation error detected: Value '{value}' at '{property_name}' failed to satisfy constraint: {specific_message}"
        super().__init__("ValidationException", message)
