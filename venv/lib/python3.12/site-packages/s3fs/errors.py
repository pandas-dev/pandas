"""S3 error codes adapted into more natural Python ones.

Adapted from: https://docs.aws.amazon.com/AmazonS3/latest/API/ErrorResponses.html
"""

import errno
import functools


# Fallback values since some systems might not have these.
ENAMETOOLONG = getattr(errno, "ENAMETOOLONG", errno.EINVAL)
ENOTEMPTY = getattr(errno, "ENOTEMPTY", errno.EINVAL)
EMSGSIZE = getattr(errno, "EMSGSIZE", errno.EINVAL)
EREMOTEIO = getattr(errno, "EREMOTEIO", errno.EIO)
EREMCHG = getattr(errno, "EREMCHG", errno.ENOENT)


ERROR_CODE_TO_EXCEPTION = {
    "AccessDenied": PermissionError,
    "AccountProblem": PermissionError,
    "AllAccessDisabled": PermissionError,
    "AmbiguousGrantByEmailAddress": functools.partial(IOError, errno.EINVAL),
    "AuthorizationHeaderMalformed": functools.partial(IOError, errno.EINVAL),
    "BadDigest": functools.partial(IOError, errno.EINVAL),
    "BucketAlreadyExists": FileExistsError,
    "BucketAlreadyOwnedByYou": FileExistsError,
    "BucketNotEmpty": functools.partial(IOError, ENOTEMPTY),
    "CredentialsNotSupported": functools.partial(IOError, errno.EINVAL),
    "CrossLocationLoggingProhibited": PermissionError,
    "EntityTooSmall": functools.partial(IOError, errno.EINVAL),
    "EntityTooLarge": functools.partial(IOError, EMSGSIZE),
    "ExpiredToken": PermissionError,
    "IllegalLocationConstraintException": PermissionError,
    "IllegalVersioningConfigurationException": functools.partial(IOError, errno.EINVAL),
    "IncompleteBody": functools.partial(IOError, errno.EINVAL),
    "IncorrectNumberOfFilesInPostRequest": functools.partial(IOError, errno.EINVAL),
    "InlineDataTooLarge": functools.partial(IOError, EMSGSIZE),
    "InternalError": functools.partial(IOError, EREMOTEIO),
    "InvalidAccessKeyId": PermissionError,
    "InvalidAddressingHeader": functools.partial(IOError, errno.EINVAL),
    "InvalidArgument": functools.partial(IOError, errno.EINVAL),
    "InvalidBucketName": functools.partial(IOError, errno.EINVAL),
    "InvalidBucketState": functools.partial(IOError, errno.EPERM),
    "InvalidDigest": functools.partial(IOError, errno.EINVAL),
    "InvalidEncryptionAlgorithmError": functools.partial(IOError, errno.EINVAL),
    "InvalidLocationConstraint": functools.partial(IOError, errno.EINVAL),
    "InvalidObjectState": PermissionError,
    "InvalidPart": functools.partial(IOError, errno.EINVAL),
    "InvalidPartOrder": functools.partial(IOError, errno.EINVAL),
    "InvalidPayer": PermissionError,
    "InvalidPolicyDocument": functools.partial(IOError, errno.EINVAL),
    "InvalidRange": functools.partial(IOError, errno.EINVAL),
    "InvalidRequest": functools.partial(IOError, errno.EINVAL),
    "InvalidSecurity": PermissionError,
    "InvalidSOAPRequest": functools.partial(IOError, errno.EINVAL),
    "InvalidStorageClass": functools.partial(IOError, errno.EINVAL),
    "InvalidTargetBucketForLogging": functools.partial(IOError, errno.EINVAL),
    "InvalidToken": functools.partial(IOError, errno.EINVAL),
    "InvalidURI": functools.partial(IOError, errno.EINVAL),
    "KeyTooLongError": functools.partial(IOError, ENAMETOOLONG),
    "MalformedACLError": functools.partial(IOError, errno.EINVAL),
    "MalformedPOSTRequest": functools.partial(IOError, errno.EINVAL),
    "MalformedXML": functools.partial(IOError, errno.EINVAL),
    "MaxMessageLengthExceeded": functools.partial(IOError, EMSGSIZE),
    "MaxPostPreDataLengthExceededError": functools.partial(IOError, EMSGSIZE),
    "MetadataTooLarge": functools.partial(IOError, EMSGSIZE),
    "MethodNotAllowed": functools.partial(IOError, errno.EPERM),
    "MissingAttachment": functools.partial(IOError, errno.EINVAL),
    "MissingContentLength": functools.partial(IOError, errno.EINVAL),
    "MissingRequestBodyError": functools.partial(IOError, errno.EINVAL),
    "MissingSecurityElement": functools.partial(IOError, errno.EINVAL),
    "MissingSecurityHeader": functools.partial(IOError, errno.EINVAL),
    "NoLoggingStatusForKey": functools.partial(IOError, errno.EINVAL),
    "NoSuchBucket": FileNotFoundError,
    "NoSuchBucketPolicy": FileNotFoundError,
    "NoSuchKey": FileNotFoundError,
    "NoSuchLifecycleConfiguration": FileNotFoundError,
    "NoSuchUpload": FileNotFoundError,
    "NoSuchVersion": FileNotFoundError,
    "NotImplemented": functools.partial(IOError, errno.ENOSYS),
    "NotSignedUp": PermissionError,
    "OperationAborted": functools.partial(IOError, errno.EBUSY),
    "PermanentRedirect": functools.partial(IOError, EREMCHG),
    "PreconditionFailed": functools.partial(IOError, errno.EINVAL),
    "Redirect": functools.partial(IOError, EREMCHG),
    "RestoreAlreadyInProgress": functools.partial(IOError, errno.EBUSY),
    "RequestIsNotMultiPartContent": functools.partial(IOError, errno.EINVAL),
    "RequestTimeout": TimeoutError,
    "RequestTimeTooSkewed": PermissionError,
    "RequestTorrentOfBucketError": functools.partial(IOError, errno.EPERM),
    "SignatureDoesNotMatch": PermissionError,
    "ServiceUnavailable": functools.partial(IOError, errno.EBUSY),
    "SlowDown": functools.partial(IOError, errno.EBUSY),
    "TemporaryRedirect": functools.partial(IOError, EREMCHG),
    "TokenRefreshRequired": functools.partial(IOError, errno.EINVAL),
    "TooManyBuckets": functools.partial(IOError, errno.EINVAL),
    "UnexpectedContent": functools.partial(IOError, errno.EINVAL),
    "UnresolvableGrantByEmailAddress": functools.partial(IOError, errno.EINVAL),
    "UserKeyMustBeSpecified": functools.partial(IOError, errno.EINVAL),
    "301": functools.partial(IOError, EREMCHG),  # PermanentRedirect
    "307": functools.partial(IOError, EREMCHG),  # Redirect
    "400": functools.partial(IOError, errno.EINVAL),
    "403": PermissionError,
    "404": FileNotFoundError,
    "405": functools.partial(IOError, errno.EPERM),
    "409": functools.partial(IOError, errno.EBUSY),
    "412": functools.partial(IOError, errno.EINVAL),  # PreconditionFailed
    "416": functools.partial(IOError, errno.EINVAL),  # InvalidRange
    "500": functools.partial(IOError, EREMOTEIO),  # InternalError
    "501": functools.partial(IOError, errno.ENOSYS),  # NotImplemented
    "503": functools.partial(IOError, errno.EBUSY),  # SlowDown
}


def translate_boto_error(error, message=None, set_cause=True, *args, **kwargs):
    """Convert a ClientError exception into a Python one.

    Parameters
    ----------

    error : botocore.exceptions.ClientError
        The exception returned by the boto API.
    message : str
        An error message to use for the returned exception. If not given, the
        error message returned by the server is used instead.
    set_cause : bool
        Whether to set the __cause__ attribute to the previous exception if the
        exception is translated.
    *args, **kwargs :
        Additional arguments to pass to the exception constructor, after the
        error message. Useful for passing the filename arguments to ``IOError``.

    Returns
    -------

    An instantiated exception ready to be thrown. If the error code isn't
    recognized, an IOError with the original error message is returned.
    """
    error_response = getattr(error, "response", None)

    if error_response is None:
        # non-http error, or response is None:
        return error
    code = error_response["Error"].get("Code")
    if (
        code == "PreconditionFailed"
        and error_response["Error"].get("Condition", "") == "If-None-Match"
    ):
        constructor = FileExistsError
    else:
        constructor = ERROR_CODE_TO_EXCEPTION.get(code)
    if constructor:
        if not message:
            message = error_response["Error"].get("Message", str(error))
        custom_exc = constructor(message, *args, **kwargs)
    else:
        # No match found, wrap this in an IOError with the appropriate message.
        custom_exc = IOError(errno.EIO, message or str(error), *args)

    if set_cause:
        custom_exc.__cause__ = error
    return custom_exc
