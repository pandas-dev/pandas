from moto.core.exceptions import ServiceException


class SNSException(ServiceException):
    pass


class SNSNotFoundError(SNSException):
    code = "NotFound"


class TopicNotFound(SNSNotFoundError):
    message = "Topic does not exist"


class ResourceNotFoundError(SNSException):
    code = "ResourceNotFound"
    message = "Resource does not exist"


class DuplicateSnsEndpointError(SNSException):
    code = "InvalidParameter"


class SnsEndpointDisabled(SNSException):
    code = "EndpointDisabled"


class SNSInvalidParameter(SNSException):
    code = "InvalidParameter"


class InvalidParameterValue(SNSException):
    code = "InvalidParameter"


class TagLimitExceededError(SNSException):
    code = "TagLimitExceeded"
    message = "Could not complete request: tag quota of per resource exceeded"


class InternalError(SNSException):
    code = "InternalError"


class TooManyEntriesInBatchRequest(SNSException):
    code = "TooManyEntriesInBatchRequest"
    message = "The batch request contains more entries than permissible."


class BatchEntryIdsNotDistinct(SNSException):
    code = "BatchEntryIdsNotDistinct"
    message = "Two or more batch entries in the request have the same Id."
