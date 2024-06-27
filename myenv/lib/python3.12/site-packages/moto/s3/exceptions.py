from typing import TYPE_CHECKING, Any, Optional, Union

from moto.core.exceptions import RESTError

from .notifications import S3NotificationEvent

if TYPE_CHECKING:
    from moto.s3.models import FakeDeleteMarker


ERROR_WITH_BUCKET_NAME = """{% extends 'single_error' %}
{% block extra %}<BucketName>{{ bucket }}</BucketName>{% endblock %}
"""

ERROR_WITH_KEY_NAME = """{% extends 'single_error' %}
{% block extra %}<Key>{{ key }}</Key>{% endblock %}
"""

ERROR_WITH_ARGUMENT = """{% extends 'single_error' %}
{% block extra %}<ArgumentName>{{ name }}</ArgumentName>
<ArgumentValue>{{ value }}</ArgumentValue>{% endblock %}
"""

ERROR_WITH_UPLOADID = """{% extends 'single_error' %}
{% block extra %}<UploadId>{{ upload_id }}</UploadId>{% endblock %}
"""

ERROR_WITH_CONDITION_NAME = """{% extends 'single_error' %}
{% block extra %}<Condition>{{ condition }}</Condition>{% endblock %}
"""

ERROR_WITH_RANGE = """{% extends 'single_error' %}
{% block extra %}<ActualObjectSize>{{ actual_size }}</ActualObjectSize>
<RangeRequested>{{ range_requested }}</RangeRequested>{% endblock %}
"""

ERROR_WITH_STORAGE_CLASS = """{% extends 'single_error' %}
{% block extra %}<StorageClass>{{ storage_class }}</StorageClass>{% endblock %}
"""


class S3ClientError(RESTError):
    # S3 API uses <RequestID> as the XML tag in response messages
    request_id_tag_name = "RequestID"

    extended_templates = {
        "bucket_error": ERROR_WITH_BUCKET_NAME,
        "key_error": ERROR_WITH_KEY_NAME,
        "argument_error": ERROR_WITH_ARGUMENT,
        "error_uploadid": ERROR_WITH_UPLOADID,
        "condition_error": ERROR_WITH_CONDITION_NAME,
        "range_error": ERROR_WITH_RANGE,
        "storage_error": ERROR_WITH_STORAGE_CLASS,
    }
    env = RESTError.extended_environment(extended_templates)

    def __init__(self, *args: Any, **kwargs: Any):
        kwargs.setdefault("template", "single_error")
        super().__init__(*args, **kwargs)


class InvalidArgumentError(S3ClientError):
    code = 400

    def __init__(self, message: str, name: str, value: str, *args: Any, **kwargs: Any):
        kwargs.setdefault("template", "argument_error")
        kwargs["name"] = name
        kwargs["value"] = value
        super().__init__("InvalidArgument", message, *args, **kwargs)


class AccessForbidden(S3ClientError):
    code = 403

    def __init__(self, msg: str):
        super().__init__("AccessForbidden", msg)


class BadRequest(S3ClientError):
    code = 403

    def __init__(self, msg: str):
        super().__init__("BadRequest", msg)


class BucketError(S3ClientError):
    def __init__(self, *args: Any, **kwargs: Any):
        kwargs.setdefault("template", "bucket_error")
        super().__init__(*args, **kwargs)


class BucketAlreadyExists(BucketError):
    code = 409

    def __init__(self, *args: Any, **kwargs: Any):
        kwargs.setdefault("template", "bucket_error")
        super().__init__(
            "BucketAlreadyExists",
            (
                "The requested bucket name is not available. The bucket "
                "namespace is shared by all users of the system. Please "
                "select a different name and try again"
            ),
            *args,
            **kwargs,
        )


class MissingBucket(BucketError):
    code = 404

    def __init__(self, bucket: str):
        super().__init__(
            "NoSuchBucket", "The specified bucket does not exist", bucket=bucket
        )


class MissingKey(S3ClientError):
    code = 404

    def __init__(self, **kwargs: Any):
        kwargs.setdefault("template", "key_error")
        super().__init__("NoSuchKey", "The specified key does not exist.", **kwargs)


class MissingVersion(S3ClientError):
    code = 404

    def __init__(self) -> None:
        super().__init__("NoSuchVersion", "The specified version does not exist.")


class InvalidVersion(S3ClientError):
    code = 400

    def __init__(self, version_id: str, *args: Any, **kwargs: Any):
        kwargs.setdefault("template", "argument_error")
        kwargs["name"] = "versionId"
        kwargs["value"] = version_id
        super().__init__(
            "InvalidArgument", "Invalid version id specified", *args, **kwargs
        )


class ObjectNotInActiveTierError(S3ClientError):
    code = 403

    def __init__(self, key_name: Any):
        super().__init__(
            "ObjectNotInActiveTierError",
            "The source object of the COPY operation is not in the active tier and is only stored in Amazon Glacier.",
            Key=key_name,
        )


class InvalidPartOrder(S3ClientError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "InvalidPartOrder",
            "The list of parts was not in ascending order. The parts list must be specified in order by part number.",
        )


class InvalidPart(S3ClientError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "InvalidPart",
            "One or more of the specified parts could not be found. The part might not have been uploaded, or the specified entity tag might not have matched the part's entity tag.",
        )


class EntityTooSmall(S3ClientError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "EntityTooSmall",
            "Your proposed upload is smaller than the minimum allowed object size.",
        )


class InvalidRequest(S3ClientError):
    code = 400

    def __init__(self, method: str):
        super().__init__(
            "InvalidRequest",
            f"Found unsupported HTTP method in CORS config. Unsupported method is {method}",
        )


class IllegalLocationConstraintException(S3ClientError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "IllegalLocationConstraintException",
            "The unspecified location constraint is incompatible for the region specific endpoint this request was sent to.",
        )


class MalformedXML(S3ClientError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "MalformedXML",
            "The XML you provided was not well-formed or did not validate against our published schema",
        )


class MalformedACLError(S3ClientError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "MalformedACLError",
            "The XML you provided was not well-formed or did not validate against our published schema",
        )


class InvalidTargetBucketForLogging(S3ClientError):
    code = 400

    def __init__(self, msg: str):
        super().__init__("InvalidTargetBucketForLogging", msg)


class CrossLocationLoggingProhibitted(S3ClientError):
    code = 403

    def __init__(self) -> None:
        super().__init__(
            "CrossLocationLoggingProhibitted", "Cross S3 location logging not allowed."
        )


class InvalidMaxPartArgument(S3ClientError):
    code = 400

    def __init__(self, arg: str, min_val: int, max_val: int):
        error = f"Argument {arg} must be an integer between {min_val} and {max_val}"
        super().__init__("InvalidArgument", error)


class InvalidMaxPartNumberArgument(InvalidArgumentError):
    code = 400

    def __init__(self, value: int):
        error = "Part number must be an integer between 1 and 10000, inclusive"
        super().__init__(message=error, name="partNumber", value=value)  # type: ignore


class NotAnIntegerException(InvalidArgumentError):
    code = 400

    def __init__(self, name: str, value: int):
        error = f"Provided {name} not an integer or within integer range"
        super().__init__(message=error, name=name, value=value)  # type: ignore


class InvalidNotificationARN(S3ClientError):
    code = 400

    def __init__(self) -> None:
        super().__init__("InvalidArgument", "The ARN is not well formed")


class InvalidNotificationDestination(S3ClientError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "InvalidArgument",
            "The notification destination service region is not valid for the bucket location constraint",
        )


class InvalidNotificationEvent(S3ClientError):
    code = 400

    def __init__(self, event_name: str) -> None:
        super().__init__(
            "InvalidArgument",
            (
                f"The event '{event_name}' is not supported for notifications. "
                f"Supported events are as follows: {S3NotificationEvent.events()}"
            ),
        )


class InvalidStorageClass(S3ClientError):
    code = 400

    def __init__(self, storage: Optional[str]):
        super().__init__(
            "InvalidStorageClass",
            "The storage class you specified is not valid",
            storage=storage,
        )


class InvalidBucketName(S3ClientError):
    code = 400

    def __init__(self) -> None:
        super().__init__("InvalidBucketName", "The specified bucket is not valid.")


class DuplicateTagKeys(S3ClientError):
    code = 400

    def __init__(self) -> None:
        super().__init__("InvalidTag", "Cannot provide multiple Tags with the same key")


class S3AccessDeniedError(S3ClientError):
    code = 403

    def __init__(self) -> None:
        super().__init__("AccessDenied", "Access Denied")


class BucketAccessDeniedError(BucketError):
    code = 403

    def __init__(self, bucket: str):
        super().__init__("AccessDenied", "Access Denied", bucket=bucket)


class S3InvalidTokenError(S3ClientError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "InvalidToken", "The provided token is malformed or otherwise invalid."
        )


class S3AclAndGrantError(S3ClientError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "InvalidRequest",
            "Specifying both Canned ACLs and Header Grants is not allowed",
        )


class BucketInvalidTokenError(BucketError):
    code = 400

    def __init__(self, bucket: str):
        super().__init__(
            "InvalidToken",
            "The provided token is malformed or otherwise invalid.",
            bucket=bucket,
        )


class S3InvalidAccessKeyIdError(S3ClientError):
    code = 403

    def __init__(self) -> None:
        super().__init__(
            "InvalidAccessKeyId",
            "The AWS Access Key Id you provided does not exist in our records.",
        )


class BucketInvalidAccessKeyIdError(S3ClientError):
    code = 403

    def __init__(self, bucket: str):
        super().__init__(
            "InvalidAccessKeyId",
            "The AWS Access Key Id you provided does not exist in our records.",
            bucket=bucket,
        )


class S3SignatureDoesNotMatchError(S3ClientError):
    code = 403

    def __init__(self) -> None:
        super().__init__(
            "SignatureDoesNotMatch",
            "The request signature we calculated does not match the signature you provided. Check your key and signing method.",
        )


class BucketSignatureDoesNotMatchError(S3ClientError):
    code = 403

    def __init__(self, bucket: str):
        super().__init__(
            "SignatureDoesNotMatch",
            "The request signature we calculated does not match the signature you provided. Check your key and signing method.",
            bucket=bucket,
        )


class NoSuchPublicAccessBlockConfiguration(S3ClientError):
    code = 404

    def __init__(self) -> None:
        super().__init__(
            "NoSuchPublicAccessBlockConfiguration",
            "The public access block configuration was not found",
        )


class InvalidPublicAccessBlockConfiguration(S3ClientError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "InvalidRequest",
            "Must specify at least one configuration.",
        )


class WrongPublicAccessBlockAccountIdError(S3ClientError):
    code = 403

    def __init__(self) -> None:
        super().__init__("AccessDenied", "Access Denied")


class NoSystemTags(S3ClientError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "InvalidTag", "System tags cannot be added/updated by requester"
        )


class NoSuchUpload(S3ClientError):
    code = 404

    def __init__(self, upload_id: Union[int, str], *args: Any, **kwargs: Any):
        kwargs.setdefault("template", "error_uploadid")
        kwargs["upload_id"] = upload_id
        super().__init__(
            "NoSuchUpload",
            "The specified upload does not exist. The upload ID may be invalid, or the upload may have been aborted or completed.",
            *args,
            **kwargs,
        )


class PreconditionFailed(S3ClientError):
    code = 412

    def __init__(self, failed_condition: str, **kwargs: Any):
        kwargs.setdefault("template", "condition_error")
        super().__init__(
            "PreconditionFailed",
            "At least one of the pre-conditions you specified did not hold",
            condition=failed_condition,
            **kwargs,
        )


class InvalidRange(S3ClientError):
    code = 416

    def __init__(self, range_requested: str, actual_size: str, **kwargs: Any):
        kwargs.setdefault("template", "range_error")
        super().__init__(
            "InvalidRange",
            "The requested range is not satisfiable",
            range_requested=range_requested,
            actual_size=actual_size,
            **kwargs,
        )


class InvalidContinuationToken(S3ClientError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "InvalidArgument", "The continuation token provided is incorrect"
        )


class InvalidObjectState(BucketError):
    code = 403

    def __init__(self, storage_class: Optional[str], **kwargs: Any):
        kwargs.setdefault("template", "storage_error")
        super().__init__(
            error_type="InvalidObjectState",
            message="The operation is not valid for the object's storage class",
            storage_class=storage_class,
            **kwargs,
        )


class LockNotEnabled(S3ClientError):
    code = 400

    def __init__(self) -> None:
        super().__init__("InvalidRequest", "Bucket is missing ObjectLockConfiguration")


class AccessDeniedByLock(S3ClientError):
    code = 400

    def __init__(self) -> None:
        super().__init__("AccessDenied", "Access Denied")


class InvalidContentMD5(S3ClientError):
    code = 400

    def __init__(self) -> None:
        super().__init__("InvalidContentMD5", "Content MD5 header is invalid")


class BucketNeedsToBeNew(S3ClientError):
    code = 400

    def __init__(self) -> None:
        super().__init__("InvalidBucket", "Bucket needs to be empty")


class BucketMustHaveLockeEnabled(S3ClientError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "InvalidBucketState",
            "Object Lock configuration cannot be enabled on existing buckets",
        )


class CopyObjectMustChangeSomething(S3ClientError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "InvalidRequest",
            "This copy request is illegal because it is trying to copy an object to itself without changing the object's metadata, storage class, website redirect location or encryption attributes.",
        )


class InvalidFilterRuleName(InvalidArgumentError):
    code = 400

    def __init__(self, value: str):
        super().__init__(
            "filter rule name must be either prefix or suffix",
            "FilterRule.Name",
            value,
        )


class InvalidTagError(S3ClientError):
    code = 400

    def __init__(self, value: str):
        super().__init__("InvalidTag", value)


class ObjectLockConfigurationNotFoundError(S3ClientError):
    code = 404

    def __init__(self) -> None:
        super().__init__(
            "ObjectLockConfigurationNotFoundError",
            "Object Lock configuration does not exist for this bucket",
        )


class HeadOnDeleteMarker(Exception):
    """Marker to indicate that we've called `head_object()` on a FakeDeleteMarker"""

    def __init__(self, marker: "FakeDeleteMarker"):
        self.marker = marker


class DaysMustNotProvidedForSelectRequest(S3ClientError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "DaysMustNotProvidedForSelectRequest",
            "`Days` must not be provided for select requests",
        )


class DaysMustProvidedExceptForSelectRequest(S3ClientError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "DaysMustProvidedExceptForSelectRequest",
            "`Days` must be provided except for select requests",
        )
