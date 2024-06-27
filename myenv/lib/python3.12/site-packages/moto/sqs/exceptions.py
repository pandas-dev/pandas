from moto.core.exceptions import RESTError


class ReceiptHandleIsInvalid(RESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "ReceiptHandleIsInvalid", "The input receipt handle is invalid."
        )


class MessageAttributesInvalid(RESTError):
    code = 400

    def __init__(self, description: str):
        super().__init__("MessageAttributesInvalid", description)


class QueueDoesNotExist(RESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "AWS.SimpleQueueService.NonExistentQueue",
            "The specified queue does not exist for this wsdl version.",
            template="wrapped_single_error",
        )


class QueueAlreadyExists(RESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("QueueAlreadyExists", message)


class EmptyBatchRequest(RESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "EmptyBatchRequest",
            "There should be at least one SendMessageBatchRequestEntry in the request.",
        )


class InvalidBatchEntryId(RESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "InvalidBatchEntryId",
            "A batch entry id can only contain alphanumeric characters, "
            "hyphens and underscores. It can be at most 80 letters long.",
        )


class BatchRequestTooLong(RESTError):
    code = 400

    def __init__(self, length: int):
        super().__init__(
            "BatchRequestTooLong",
            "Batch requests cannot be longer than 262144 bytes. "
            f"You have sent {length} bytes.",
        )


class BatchEntryIdsNotDistinct(RESTError):
    code = 400

    def __init__(self, entry_id: str):
        super().__init__("BatchEntryIdsNotDistinct", f"Id {entry_id} repeated.")


class TooManyEntriesInBatchRequest(RESTError):
    code = 400

    def __init__(self, number: int):
        super().__init__(
            "TooManyEntriesInBatchRequest",
            "Maximum number of entries per request are 10. " f"You have sent {number}.",
        )


class InvalidAttributeName(RESTError):
    code = 400

    def __init__(self, attribute_name: str):
        super().__init__("InvalidAttributeName", f"Unknown Attribute {attribute_name}.")


class InvalidAttributeValue(RESTError):
    code = 400

    def __init__(self, attribute_name: str):
        super().__init__(
            "InvalidAttributeValue",
            f"Invalid value for the parameter {attribute_name}.",
        )


class InvalidParameterValue(RESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("InvalidParameterValue", message)


class MissingParameter(RESTError):
    code = 400

    def __init__(self, parameter: str):
        super().__init__(
            "MissingParameter", f"The request must contain the parameter {parameter}."
        )


class OverLimit(RESTError):
    code = 403

    def __init__(self, count: int):
        super().__init__(
            "OverLimit", f"{count} Actions were found, maximum allowed is 7."
        )
