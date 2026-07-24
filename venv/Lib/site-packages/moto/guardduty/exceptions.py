from typing import Any

from moto.core.exceptions import JsonRESTError


class GuardDutyException(JsonRESTError):
    pass


class DetectorNotFoundException(GuardDutyException):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "InvalidInputException",
            "The request is rejected because the input detectorId is not owned by the current account.",
        )

    def get_headers(self, *args: Any, **kwargs: Any) -> list[tuple[str, str]]:
        return [("X-Amzn-ErrorType", "BadRequestException")]


class FilterNotFoundException(GuardDutyException):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "InvalidInputException",
            "The request is rejected since no such resource found.",
        )

    def get_headers(self, *args: Any, **kwargs: Any) -> list[tuple[str, str]]:
        return [("X-Amzn-ErrorType", "BadRequestException")]
