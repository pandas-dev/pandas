""" Note: this code is used only by the static type checker!

_see: <https://gist.github.com/alexcasalboni/a545b68ee164b165a74a20a5fee9d133>
_and: <https://docs.aws.amazon.com/lambda/latest/dg/python-context-object.html>

"""
from typing import Dict, Any

LambdaDict = Dict[str, Any]


class LambdaCognitoIdentity(object):
    cognito_identity_id: str
    cognito_identity_pool_id: str


class LambdaClientContextMobileClient(object):
    installation_id: str
    app_title: str
    app_version_name: str
    app_version_code: str
    app_package_name: str


class LambdaClientContext(object):
    client: LambdaClientContextMobileClient
    custom: LambdaDict
    env: LambdaDict


class LambdaContext(object):
    function_name: str
    function_version: str
    invoked_function_arn: str
    memory_limit_in_mb: int
    aws_request_id: str
    log_group_name: str
    log_stream_name: str
    deadline_ms: int
    identity: LambdaCognitoIdentity
    client_context: LambdaClientContext

    @staticmethod
    def get_remaining_time_in_millis() -> int:
        return 0
