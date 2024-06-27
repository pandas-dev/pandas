import json
from json import JSONDecodeError
from typing import Any, Final, Optional

from moto.stepfunctions.parser.asl.eval.environment import Environment
from moto.stepfunctions.parser.asl.utils.boto_client import boto_client_for
from moto.stepfunctions.parser.asl.utils.encoding import to_json_str


class LambdaFunctionErrorException(Exception):
    function_error: Final[Optional[str]]
    payload: Final[str]

    def __init__(self, function_error: Optional[str], payload: str):
        self.function_error = function_error
        self.payload = payload


def exec_lambda_function(
    env: Environment, parameters: dict, region: str, account: str
) -> None:
    lambda_client = boto_client_for(region=region, account=account, service="lambda")

    invocation_resp = lambda_client.invoke(**parameters)

    func_error: Optional[str] = invocation_resp.get("FunctionError")
    payload_json = json.load(invocation_resp["Payload"])
    if func_error:
        payload_str = json.dumps(payload_json, separators=(",", ":"))
        raise LambdaFunctionErrorException(func_error, payload_str)

    invocation_resp["Payload"] = payload_json

    env.stack.append(invocation_resp)


def to_payload_type(payload: Any) -> Optional[bytes]:
    if isinstance(payload, bytes):
        return payload

    if payload is None:
        str_value = to_json_str(dict())
    elif isinstance(payload, str):
        try:
            json.loads(payload)
            str_value = payload
        except JSONDecodeError:
            str_value = to_json_str(payload)
    else:
        str_value = to_json_str(payload)
    return str_value.encode("utf-8")
