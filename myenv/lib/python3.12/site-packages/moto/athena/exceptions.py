import json

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
