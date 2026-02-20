import json
from typing import Optional

from moto.core.exceptions import JsonRESTError


class AthenaClientError(JsonRESTError):
    def __init__(self, code: str, message: str):
        super().__init__(error_type="InvalidRequestException", message=message)
        self.description = json.dumps(
            {
                "Error": {
                    "Code": code,
                    "Message": message,
                    "Type": "InvalidRequestException",
                },
                "RequestId": "6876f774-7273-11e4-85dc-39e55ca848d1",
            }
        )


class InvalidArgumentException(JsonRESTError):
    """The specified input parameter has a value that is not valid."""

    code = 400

    def __init__(self, message: str):
        super().__init__("InvalidArgumentException", message)


class QueryStillRunning(JsonRESTError):
    def __init__(self, current_status: Optional[str]):
        msg = f"Query has not yet finished. Current state: {current_status}"
        super().__init__("InvalidRequestException", msg)
