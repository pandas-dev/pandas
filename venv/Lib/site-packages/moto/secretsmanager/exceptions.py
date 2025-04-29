from moto.core.exceptions import JsonRESTError


class SecretsManagerClientError(JsonRESTError):
    code = 400


class ResourceNotFoundException(SecretsManagerClientError):
    def __init__(self, message: str):
        self.code = 404
        super().__init__("ResourceNotFoundException", message)


class SecretNotFoundException(SecretsManagerClientError):
    def __init__(self) -> None:
        self.code = 404
        super().__init__(
            "ResourceNotFoundException",
            message="Secrets Manager can't find the specified secret.",
        )


class SecretHasNoValueException(SecretsManagerClientError):
    def __init__(self, version_stage: str):
        self.code = 404
        super().__init__(
            "ResourceNotFoundException",
            message="Secrets Manager can't find the specified secret "
            f"value for staging label: {version_stage}",
        )


class SecretStageVersionMismatchException(SecretsManagerClientError):
    def __init__(self) -> None:
        self.code = 404
        super().__init__(
            "InvalidRequestException",
            message="You provided a VersionStage that is not associated to the provided VersionId.",
        )


class ClientError(SecretsManagerClientError):
    def __init__(self, message: str):
        super().__init__("InvalidParameterValue", message)


class InvalidParameterException(SecretsManagerClientError):
    def __init__(self, message: str):
        super().__init__("InvalidParameterException", message)


class ResourceExistsException(SecretsManagerClientError):
    def __init__(self, message: str):
        super().__init__("ResourceExistsException", message)


class InvalidRequestException(SecretsManagerClientError):
    def __init__(self, message: str):
        super().__init__("InvalidRequestException", message)


class ValidationException(SecretsManagerClientError):
    def __init__(self, message: str):
        super().__init__("ValidationException", message)


class OperationNotPermittedOnReplica(InvalidParameterException):
    def __init__(self) -> None:
        super().__init__(
            "Operation not permitted on a replica secret. Call must be made in primary secret's region."
        )
