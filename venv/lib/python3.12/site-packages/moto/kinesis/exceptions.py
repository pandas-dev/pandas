from typing import Optional

from moto.core.exceptions import JsonRESTError


class ResourceNotFoundError(JsonRESTError):
    def __init__(self, message: str):
        super().__init__(error_type="ResourceNotFoundException", message=message)


class ResourceInUseError(JsonRESTError):
    def __init__(self, message: str):
        super().__init__(error_type="ResourceInUseException", message=message)


class StreamNotFoundError(ResourceNotFoundError):
    def __init__(self, stream_name: str, account_id: str):
        super().__init__(f"Stream {stream_name} under account {account_id} not found.")


class StreamCannotBeUpdatedError(JsonRESTError):
    def __init__(self, stream_name: str, account_id: str):
        message = f"Request is invalid. Stream {stream_name} under account {account_id} is in On-Demand mode."
        super().__init__(error_type="ValidationException", message=message)


class ShardNotFoundError(ResourceNotFoundError):
    def __init__(self, shard_id: str, stream: str, account_id: str):
        super().__init__(
            f"Could not find shard {shard_id} in stream {stream} under account {account_id}."
        )


class ConsumerNotFound(ResourceNotFoundError):
    def __init__(self, consumer: str, account_id: str):
        super().__init__(f"Consumer {consumer}, account {account_id} not found.")


class InvalidArgumentError(JsonRESTError):
    def __init__(self, message: str):
        super().__init__(error_type="InvalidArgumentException", message=message)


class InvalidRetentionPeriod(InvalidArgumentError):
    def __init__(self, hours: int, too_short: bool):
        if too_short:
            msg = f"Minimum allowed retention period is 24 hours. Requested retention period ({hours} hours) is too short."
        else:
            msg = f"Maximum allowed retention period is 8760 hours. Requested retention period ({hours} hours) is too long."
        super().__init__(msg)


class InvalidDecreaseRetention(InvalidArgumentError):
    def __init__(self, name: Optional[str], requested: int, existing: int):
        msg = f"Requested retention period ({requested} hours) for stream {name} can not be longer than existing retention period ({existing} hours). Use IncreaseRetentionPeriod API."
        super().__init__(msg)


class InvalidIncreaseRetention(InvalidArgumentError):
    def __init__(self, name: Optional[str], requested: int, existing: int):
        msg = f"Requested retention period ({requested} hours) for stream {name} can not be shorter than existing retention period ({existing} hours). Use DecreaseRetentionPeriod API."
        super().__init__(msg)


class ValidationException(JsonRESTError):
    def __init__(self, value: str, position: str, regex_to_match: str):
        msg = f"1 validation error detected: Value '{value}' at '{position}' failed to satisfy constraint: Member must satisfy regular expression pattern: {regex_to_match}"
        super().__init__(error_type="ValidationException", message=msg)


class RecordSizeExceedsLimit(JsonRESTError):
    def __init__(self, position: int):
        msg = f"1 validation error detected: Value at 'records.{position}.member.data' failed to satisfy constraint: Member must have length less than or equal to 1048576"
        super().__init__(error_type="ValidationException", message=msg)


class TotalRecordsSizeExceedsLimit(JsonRESTError):
    def __init__(self) -> None:
        super().__init__(
            error_type="InvalidArgumentException",
            message="Records size exceeds 5 MB limit",
        )


class TooManyRecords(JsonRESTError):
    def __init__(self) -> None:
        msg = "1 validation error detected: Value at 'records' failed to satisfy constraint: Member must have length less than or equal to 500"
        super().__init__(error_type="ValidationException", message=msg)
