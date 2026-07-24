from moto.core.exceptions import ServiceException


class CloudFrontException(ServiceException):
    pass


class OriginDoesNotExist(CloudFrontException):
    def __init__(self) -> None:
        super().__init__(
            "NoSuchOrigin",
            "One or more of your origins or origin groups do not exist.",
        )


class DomainNameNotAnS3Bucket(CloudFrontException):
    def __init__(self) -> None:
        super().__init__(
            "InvalidArgument",
            "The parameter Origin DomainName does not refer to a valid S3 bucket.",
        )


class DistributionAlreadyExists(CloudFrontException):
    def __init__(self, dist_id: str):
        super().__init__(
            "DistributionAlreadyExists",
            f"The caller reference that you are using to create a distribution is associated with another distribution. Already exists: {dist_id}",
        )


class InvalidIfMatchVersion(CloudFrontException):
    def __init__(self) -> None:
        super().__init__(
            "InvalidIfMatchVersion",
            "The If-Match version is missing or not valid for the resource.",
        )


class NoSuchDistribution(CloudFrontException):
    def __init__(self) -> None:
        super().__init__(
            "NoSuchDistribution", "The specified distribution does not exist."
        )


class NoSuchOriginAccessControl(CloudFrontException):
    def __init__(self) -> None:
        super().__init__(
            "NoSuchOriginAccessControl",
            "The specified origin access control does not exist.",
        )


class NoSuchInvalidation(CloudFrontException):
    def __init__(self) -> None:
        super().__init__(
            "NoSuchInvalidation", "The specified invalidation does not exist."
        )
