from __future__ import annotations

from typing import TYPE_CHECKING, Any

from werkzeug.exceptions import RequestedRangeNotSatisfiable

from moto.core.exceptions import ServiceException

from .notifications import S3NotificationEvent

if TYPE_CHECKING:
    from moto.s3.models import FakeDeleteMarker


class S3ClientError(ServiceException):
    pass


class InvalidArgumentError(S3ClientError):
    code = "InvalidArgument"

    def __init__(self, message: str, name: str, value: str):
        super().__init__("InvalidArgument", message)
        self.argument_name = name
        self.argument_value = value


class AccessForbidden(S3ClientError):
    def __init__(self, msg: str):
        super().__init__("AccessForbidden", msg)


class BadRequest(S3ClientError):
    def __init__(self, msg: str):
        super().__init__("BadRequest", msg)


class BucketError(S3ClientError):
    pass


class BucketAlreadyExists(BucketError):
    code = "BucketAlreadyExists"

    def __init__(self, bucket: str):
        super().__init__(
            "The requested bucket name is not available. The bucket "
            "namespace is shared by all users of the system. Please "
            "select a different name and try again"
        )
        self.bucket_name = bucket


class MissingBucket(BucketError):
    code = "NoSuchBucket"

    def __init__(self, bucket: str):
        super().__init__("The specified bucket does not exist")
        self.bucket_name = bucket


class MissingKey(S3ClientError):
    code = "NoSuchKey"

    def __init__(self, key: str | None = None):
        super().__init__("The specified key does not exist.")
        self.key = key


class MissingVersion(S3ClientError):
    code = "NoSuchVersion"

    def __init__(self, key: str | None = None, version_id: str | None = None) -> None:
        super().__init__("The specified version does not exist.")
        self.key = key
        self.version_id = version_id


class MissingInventoryConfig(S3ClientError):
    code = "NoSuchInventoryConfig"

    def __init__(self) -> None:
        super().__init__("The specified inventory configuration does not exist.")


class InvalidVersion(S3ClientError):
    code = "InvalidArgument"

    def __init__(self, version_id: str):
        super().__init__("Invalid version id specified")
        self.argument_name = "versionId"
        self.argument_value = version_id


class ObjectNotInActiveTierError(S3ClientError):
    code = "ObjectNotInActiveTierError"

    def __init__(self, key_name: Any):
        super().__init__(
            "The source object of the COPY operation is not in the active tier and is only stored in Amazon Glacier."
        )
        self.key = key_name


class InvalidPartOrder(S3ClientError):
    code = "InvalidPartOrder"

    def __init__(self) -> None:
        super().__init__(
            "The list of parts was not in ascending order. The parts list must be specified in order by part number."
        )


class InvalidPart(S3ClientError):
    code = "InvalidPart"

    def __init__(self) -> None:
        super().__init__(
            "One or more of the specified parts could not be found. The part might not have been uploaded, or the specified entity tag might not have matched the part's entity tag."
        )


class EntityTooSmall(S3ClientError):
    def __init__(self) -> None:
        super().__init__(
            "EntityTooSmall",
            "Your proposed upload is smaller than the minimum allowed object size.",
        )


class InvalidRequest(S3ClientError):
    def __init__(self, method: str):
        super().__init__(
            "InvalidRequest",
            f"Found unsupported HTTP method in CORS config. Unsupported method is {method}",
        )


class IllegalLocationConstraintException(S3ClientError):
    def __init__(self) -> None:
        super().__init__(
            "IllegalLocationConstraintException",
            "The unspecified location constraint is incompatible for the region specific endpoint this request was sent to.",
        )


class IncompatibleLocationConstraintException(S3ClientError):
    def __init__(self, location: str) -> None:
        super().__init__(
            "IllegalLocationConstraintException",
            f"The {location} location constraint is incompatible for the region specific endpoint this request was sent to.",
        )


class InvalidLocationConstraintException(S3ClientError):
    def __init__(self) -> None:
        super().__init__(
            "InvalidLocationConstraint",
            "The specified location-constraint is not valid",
        )


class InvalidNamespaceHeaderException(S3ClientError):
    def __init__(self, account_id: str, region: str) -> None:
        super().__init__(
            "InvalidNamespaceHeader",
            "The requested bucket name did not include the account-regional namespace suffix, "
            "but the provided x-amz-bucket-namespace header value is account-regional. "
            f"Specify -{account_id}-{region}-an as the bucket name suffix to create a bucket "
            "in your account-regional namespace, or remove the header.",
        )


class MalformedXML(S3ClientError):
    def __init__(self) -> None:
        super().__init__(
            "MalformedXML",
            "The XML you provided was not well-formed or did not validate against our published schema",
        )


class MalformedACLError(S3ClientError):
    def __init__(self) -> None:
        super().__init__(
            "MalformedACLError",
            "The XML you provided was not well-formed or did not validate against our published schema",
        )


class InvalidTargetBucketForLogging(S3ClientError):
    def __init__(self, msg: str):
        super().__init__("InvalidTargetBucketForLogging", msg)


class CrossLocationLoggingProhibited(S3ClientError):
    def __init__(self) -> None:
        super().__init__(
            "CrossLocationLoggingProhibited", "Cross S3 location logging not allowed."
        )


class InvalidMaxPartArgument(S3ClientError):
    def __init__(self, arg: str, min_val: int, max_val: int):
        error = f"Argument {arg} must be an integer between {min_val} and {max_val}"
        super().__init__("InvalidArgument", error)


class InvalidMaxPartNumberArgument(InvalidArgumentError):
    def __init__(self, value: int):
        error = "Part number must be an integer between 1 and 10000, inclusive"
        super().__init__(message=error, name="partNumber", value=value)  # type: ignore


class NotAnIntegerException(InvalidArgumentError):
    def __init__(self, name: str, value: int):
        error = f"Provided {name} not an integer or within integer range"
        super().__init__(message=error, name=name, value=value)  # type: ignore


class InvalidNotificationARN(S3ClientError):
    def __init__(self) -> None:
        super().__init__("InvalidArgument", "The ARN is not well formed")


class InvalidNotificationDestination(S3ClientError):
    def __init__(self) -> None:
        super().__init__(
            "InvalidArgument",
            "The notification destination service region is not valid for the bucket location constraint",
        )


class InvalidNotificationEvent(S3ClientError):
    def __init__(self, event_name: str) -> None:
        super().__init__(
            "InvalidArgument",
            (
                f"The event '{event_name}' is not supported for notifications. "
                f"Supported events are as follows: {S3NotificationEvent.events()}"
            ),
        )


class InvalidStorageClass(S3ClientError):
    def __init__(self, storage: str | None):
        super().__init__(
            "InvalidStorageClass", "The storage class you specified is not valid"
        )
        self.storage_class = storage


class InvalidBucketName(S3ClientError):
    def __init__(self) -> None:
        super().__init__("InvalidBucketName", "The specified bucket is not valid.")


class DuplicateTagKeys(S3ClientError):
    def __init__(self) -> None:
        super().__init__("InvalidTag", "Cannot provide multiple Tags with the same key")


class S3AccessDeniedError(S3ClientError):
    def __init__(self) -> None:
        super().__init__("AccessDenied", "Access Denied")


class BucketAccessDeniedError(S3AccessDeniedError):
    def __init__(self, bucket: str):
        super().__init__()
        self.bucket_name = bucket


class S3InvalidTokenError(S3ClientError):
    def __init__(self) -> None:
        super().__init__(
            "InvalidToken", "The provided token is malformed or otherwise invalid."
        )


class S3AclAndGrantError(S3ClientError):
    def __init__(self) -> None:
        super().__init__(
            "InvalidRequest",
            "Specifying both Canned ACLs and Header Grants is not allowed",
        )


class BucketInvalidTokenError(S3InvalidTokenError):
    def __init__(self, bucket: str):
        super().__init__()
        self.bucket_name = bucket


class S3InvalidAccessKeyIdError(S3ClientError):
    def __init__(self) -> None:
        super().__init__(
            "InvalidAccessKeyId",
            "The AWS Access Key Id you provided does not exist in our records.",
        )


class BucketInvalidAccessKeyIdError(S3InvalidAccessKeyIdError):
    def __init__(self, bucket: str):
        super().__init__()
        self.bucket_name = bucket


class S3SignatureDoesNotMatchError(S3ClientError):
    def __init__(self) -> None:
        super().__init__(
            "SignatureDoesNotMatch",
            "The request signature we calculated does not match the signature you provided. Check your key and signing method.",
        )


class BucketSignatureDoesNotMatchError(S3SignatureDoesNotMatchError):
    def __init__(self, bucket: str):
        super().__init__()
        self.bucket_name = bucket


class NoSuchPublicAccessBlockConfiguration(S3ClientError):
    def __init__(self) -> None:
        super().__init__(
            "NoSuchPublicAccessBlockConfiguration",
            "The public access block configuration was not found",
        )


class InvalidPublicAccessBlockConfiguration(S3ClientError):
    def __init__(self) -> None:
        super().__init__(
            "InvalidRequest",
            "Must specify at least one configuration.",
        )


class WrongPublicAccessBlockAccountIdError(S3ClientError):
    def __init__(self) -> None:
        super().__init__("AccessDenied", "Access Denied")


class NoSystemTags(S3ClientError):
    def __init__(self) -> None:
        super().__init__(
            "InvalidTag", "System tags cannot be added/updated by requester"
        )


class NoSuchUpload(S3ClientError):
    code = "NoSuchUpload"

    def __init__(self, upload_id: int | str):
        super().__init__(
            "The specified upload does not exist. The upload ID may be invalid, or the upload may have been aborted or completed."
        )
        self.upload_id = upload_id


class PreconditionFailed(S3ClientError):
    code = "PreconditionFailed"

    def __init__(self, failed_condition: str):
        super().__init__(
            "At least one of the pre-conditions you specified did not hold"
        )
        self.condition = failed_condition


class InvalidRange(S3ClientError):
    code = "InvalidRange"

    def __init__(self, range_requested: str, actual_size: str):
        super().__init__("The requested range is not satisfiable")
        self.range_requested = range_requested
        self.actual_object_size = actual_size


class RangeNotSatisfiable(RequestedRangeNotSatisfiable):
    def __init__(self) -> None:
        super().__init__(description="Requested Range Not Satisfiable")


class InvalidContinuationToken(S3ClientError):
    def __init__(self) -> None:
        super().__init__(
            "InvalidArgument", "The continuation token provided is incorrect"
        )


class InvalidBucketState(S3ClientError):
    def __init__(self, msg: str):
        super().__init__("InvalidBucketState", msg)


class InvalidObjectState(S3ClientError):
    code = "InvalidObjectState"

    def __init__(self, storage_class: str | None):
        super().__init__("The operation is not valid for the object's storage class")
        self.storage_class = storage_class


class LockNotEnabled(S3ClientError):
    def __init__(self) -> None:
        super().__init__("InvalidRequest", "Bucket is missing ObjectLockConfiguration")


class MissingRequestBody(S3ClientError):
    def __init__(self) -> None:
        super().__init__("MissingRequestBodyError", "Request Body is empty")


class AccessDeniedByLock(S3ClientError):
    def __init__(self) -> None:
        super().__init__("AccessDenied", "Access Denied")


class MissingUploadObjectWithObjectLockHeaders(S3ClientError):
    def __init__(self) -> None:
        super().__init__(
            "MissingUploadObjectWithObjectLockHeaders",
            "Content-MD5 or x-amz-sdk-checksum-algorithm header required to upload an object with a retention period configured using Object Lock",
        )


class BucketNeedsToBeNew(S3ClientError):
    def __init__(self) -> None:
        super().__init__("InvalidBucket", "Bucket needs to be empty")


class CopyObjectMustChangeSomething(S3ClientError):
    def __init__(self) -> None:
        super().__init__(
            "InvalidRequest",
            "This copy request is illegal because it is trying to copy an object to itself without changing the object's metadata, storage class, website redirect location or encryption attributes.",
        )


class InvalidFilterRuleName(InvalidArgumentError):
    def __init__(self, value: str):
        super().__init__(
            "filter rule name must be either prefix or suffix",
            "FilterRule.Name",
            value,
        )


class InvalidTagError(S3ClientError):
    def __init__(self, value: str):
        super().__init__("InvalidTag", value)


class ObjectLockConfigurationNotFoundError(S3ClientError):
    def __init__(self) -> None:
        super().__init__(
            "ObjectLockConfigurationNotFoundError",
            "Object Lock configuration does not exist for this bucket",
        )


class HeadOnDeleteMarker(Exception):
    """Marker to indicate that we've called `head_object()` on a FakeDeleteMarker"""

    def __init__(self, marker: FakeDeleteMarker):
        self.marker = marker


class MethodNotAllowed(S3ClientError):
    code = "MethodNotAllowed"

    def __init__(self, method: str | None = None, resource_type: str | None = None):
        super().__init__("The specified method is not allowed against this resource.")
        self.method = method
        self.resource_type = resource_type


class NoSuchLifecycleConfiguration(S3ClientError):
    code = "NoSuchLifecycleConfiguration"

    def __init__(self, bucket_name: str):
        super().__init__("The lifecycle configuration does not exist")
        self.bucket_name = bucket_name


class NoSuchCORSConfiguration(S3ClientError):
    code = "NoSuchCORSConfiguration"

    def __init__(self, bucket_name: str):
        super().__init__("The CORS configuration does not exist")
        self.bucket_name = bucket_name


class OwnershipControlsNotFoundError(S3ClientError):
    code = "OwnershipControlsNotFoundError"

    def __init__(self, bucket_name: str):
        super().__init__("The bucket ownership controls were not found")
        self.bucket_name = bucket_name


class BucketNotEmpty(S3ClientError):
    code = "BucketNotEmpty"

    def __init__(self, bucket_name: str):
        super().__init__("The bucket you tried to delete is not empty")
        self.bucket_name = bucket_name


class NoSuchBucketPolicy(S3ClientError):
    code = "NoSuchBucketPolicy"

    def __init__(self, bucket_name: str):
        super().__init__("The bucket policy does not exist")
        self.bucket_name = bucket_name


class NoSuchTagSet(S3ClientError):
    code = "NoSuchTagSet"

    def __init__(self, bucket_name: str):
        super().__init__("The TagSet does not exist")
        self.bucket_name = bucket_name


class NoSuchWebsiteConfiguration(S3ClientError):
    code = "NoSuchWebsiteConfiguration"

    def __init__(self, bucket_name: str):
        super().__init__("The specified bucket does not have a website configuration")
        self.bucket_name = bucket_name


class ServerSideEncryptionConfigurationNotFoundError(S3ClientError):
    code = "ServerSideEncryptionConfigurationNotFoundError"

    def __init__(self, bucket_name: str):
        super().__init__("The server side encryption configuration was not found")
        self.bucket_name = bucket_name


class BucketAlreadyOwnedByYou(S3ClientError):
    code = "BucketAlreadyOwnedByYou"

    def __init__(self, bucket_name: str):
        super().__init__(
            "Your previous request to create the named bucket succeeded and you already own it."
        )
        self.bucket_name = bucket_name


class ReplicationConfigurationNotFoundError(S3ClientError):
    code = "ReplicationConfigurationNotFoundError"

    def __init__(self, bucket_name: str):
        super().__init__("The replication configuration was not found")
        self.bucket_name = bucket_name


class VersioningNotEnabledForReplication(S3ClientError):
    code = "InvalidRequest"

    def __init__(self, bucket_name: str):
        super().__init__(
            "Versioning must be 'Enabled' on the bucket to apply a replication configuration"
        )
        self.bucket_name = bucket_name
