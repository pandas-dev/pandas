import json
from typing import Optional

from moto.core.exceptions import JsonRESTError


class IoTClientError(JsonRESTError):
    code = 400


class ResourceNotFoundException(IoTClientError):
    def __init__(self, msg: Optional[str] = None):
        self.code = 404
        super().__init__(
            "ResourceNotFoundException", msg or "The specified resource does not exist"
        )


class InvalidRequestException(IoTClientError):
    def __init__(self, msg: Optional[str] = None):
        self.code = 400
        super().__init__("InvalidRequestException", msg or "The request is not valid.")


class InvalidStateTransitionException(IoTClientError):
    def __init__(self, msg: Optional[str] = None):
        self.code = 409
        super().__init__(
            "InvalidStateTransitionException",
            msg or "An attempt was made to change to an invalid state.",
        )


class VersionConflictException(IoTClientError):
    def __init__(self, name: str):
        self.code = 409
        super().__init__(
            "VersionConflictException",
            f"The version for thing {name} does not match the expected version.",
        )


class CertificateStateException(IoTClientError):
    def __init__(self, msg: str, cert_id: str):
        self.code = 406
        super().__init__("CertificateStateException", f"{msg} Id: {cert_id}")


class DeleteConflictException(IoTClientError):
    def __init__(self, msg: str):
        self.code = 409
        super().__init__("DeleteConflictException", msg)


class ResourceAlreadyExistsException(IoTClientError):
    def __init__(self, msg: str, resource_id: str, resource_arn: str):
        self.code = 409
        super().__init__(
            "ResourceAlreadyExistsException", msg or "The resource already exists."
        )
        self.description = json.dumps(
            {
                "message": self.message,
                "resourceId": resource_id,
                "resourceArn": resource_arn,
            }
        )


class VersionsLimitExceededException(IoTClientError):
    def __init__(self, name: str):
        self.code = 409
        super().__init__(
            "VersionsLimitExceededException",
            f"The policy {name} already has the maximum number of versions (5)",
        )


class ThingStillAttached(IoTClientError):
    def __init__(self, name: str):
        self.code = 409
        super().__init__(
            "InvalidRequestException",
            f"Cannot delete. Thing {name} is still attached to one or more principals",
        )
