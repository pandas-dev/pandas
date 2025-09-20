import re
from typing import List, Optional, TypedDict


class _docker_config(TypedDict, total=False):
    use_docker: bool


class _passthrough_config(TypedDict, total=False):
    services: List[str]
    urls: List[str]


class _core_config(TypedDict, total=False):
    mock_credentials: bool
    passthrough: _passthrough_config
    reset_boto3_session: bool
    service_whitelist: Optional[List[str]]


class _iam_config(TypedDict, total=False):
    load_aws_managed_policies: bool


class _sfn_config(TypedDict, total=False):
    execute_state_machine: bool


class _iot_config(TypedDict, total=False):
    use_valid_cert: bool


DefaultConfig = TypedDict(
    "DefaultConfig",
    {
        "batch": _docker_config,
        "core": _core_config,
        "lambda": _docker_config,
        "iam": _iam_config,
        "stepfunctions": _sfn_config,
        "iot": _iot_config,
    },
    total=False,
)

default_user_config: DefaultConfig = {
    "batch": {"use_docker": True},
    "lambda": {"use_docker": True},
    "core": {
        "mock_credentials": True,
        "passthrough": {"urls": [], "services": []},
        "reset_boto3_session": True,
        "service_whitelist": None,
    },
    "iam": {"load_aws_managed_policies": False},
    "stepfunctions": {"execute_state_machine": False},
    "iot": {"use_valid_cert": False},
}


def service_whitelisted(service: str) -> bool:
    services_whitelisted = default_user_config.get("core", {}).get("service_whitelist")
    return services_whitelisted is None or service in services_whitelisted


def passthrough_service(service: str) -> bool:
    passthrough_services = (
        default_user_config.get("core", {}).get("passthrough", {}).get("services", [])
    )
    return service in passthrough_services


def passthrough_url(clean_url: str) -> bool:
    passthrough_urls = (
        default_user_config.get("core", {}).get("passthrough", {}).get("urls", [])
    )
    return any([re.match(url, clean_url) for url in passthrough_urls])


def mock_credentials() -> bool:
    return (
        default_user_config.get("core", {}).get("mock_credentials", True) is not False
    )
