from moto.core.exceptions import AWSError


class AWSValidationException(AWSError):
    TYPE = "ValidationException"


class AWSResourceNotFoundException(AWSError):
    TYPE = "ResourceNotFoundException"


class CertificateNotFound(AWSResourceNotFoundException):
    def __init__(self, arn: str, account_id: str):
        super().__init__(
            message=f"Certificate with arn {arn} not found in account {account_id}"
        )


class AWSTooManyTagsException(AWSError):
    TYPE = "TooManyTagsException"
