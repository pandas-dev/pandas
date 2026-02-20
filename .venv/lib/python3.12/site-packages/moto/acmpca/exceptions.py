"""Exceptions raised by the acmpca service."""

from moto.core.exceptions import JsonRESTError


class ResourceNotFoundException(JsonRESTError):
    def __init__(self, arn: str):
        super().__init__("ResourceNotFoundException", f"Resource {arn} not found")


class InvalidS3ObjectAclInCrlConfiguration(JsonRESTError):
    code = 400

    def __init__(self, value: str):
        super().__init__(
            "InvalidS3ObjectAclInCrlConfiguration",
            f"Invalid value for parameter RevocationConfiguration.CrlConfiguration.S3ObjectAcl, value: {value}, valid values: ['PUBLIC_READ', 'BUCKET_OWNER_FULL_CONTROL']",
        )


class InvalidStateException(JsonRESTError):
    code = 400

    def __init__(self, arn: str):
        super().__init__(
            "InvalidStateException",
            f"The certificate authority {arn} is not in the correct state to have a certificate signing request.",
        )


class MalformedCertificateAuthorityException(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "MalformedCertificateAuthorityException",
            "Malformed certificate.",
        )


class InvalidPolicyException(JsonRESTError):
    def __init__(self) -> None:
        super().__init__(
            "InvalidPolicyException",
            "The resource policy is invalid or is missing a required statement.",
        )


class LockoutPreventedException(JsonRESTError):
    def __init__(self) -> None:
        super().__init__(
            "LockoutPreventedException",
            "The current action was prevented because it would lock the caller out from performing subsequent actions.",
        )


class ConcurrentModificationException(JsonRESTError):
    def __init__(self) -> None:
        super().__init__(
            "ConcurrentModificationException",
            "A previous update to your private CA is still ongoing.",
        )


class RequestInProgressException(JsonRESTError):
    def __init__(self, message: str):
        super().__init__(
            "RequestInProgressException",
            message,
        )
