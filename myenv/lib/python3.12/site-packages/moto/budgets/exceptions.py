"""Exceptions raised by the budgets service."""

from moto.core.exceptions import JsonRESTError


class DuplicateRecordException(JsonRESTError):
    code = 400

    def __init__(self, record_type: str, record_name: str):
        super().__init__(
            __class__.__name__,  # type: ignore[name-defined]
            f"Error creating {record_type}: {record_name} - the {record_type} already exists.",
        )


class NotFoundException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__(__class__.__name__, message)  # type: ignore[name-defined]


class BudgetMissingLimit(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "InvalidParameterException",
            "Unable to create/update budget - please provide one of the followings: Budget Limit/ Planned Budget Limit/ Auto Adjust Data",
        )
