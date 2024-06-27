from collections import namedtuple
from functools import partial
from typing import TYPE_CHECKING, Any, Callable

from moto.utilities.utils import PARTITION_NAMES, get_partition

if TYPE_CHECKING:
    from .models import LambdaBackend

ARN = namedtuple("ARN", ["region", "account", "function_name", "version"])
LAYER_ARN = namedtuple("LAYER_ARN", ["region", "account", "layer_name", "version"])


def make_arn(resource_type: str, region: str, account: str, name: str) -> str:
    return (
        f"arn:{get_partition(region)}:lambda:{region}:{account}:{resource_type}:{name}"
    )


make_function_arn = partial(make_arn, "function")
make_layer_arn = partial(make_arn, "layer")


def make_ver_arn(
    resource_type: str, region: str, account: str, name: str, version: str = "1"
) -> str:
    arn = make_arn(resource_type, region, account, name)
    return f"{arn}:{version}"


make_function_ver_arn = partial(make_ver_arn, "function")
make_layer_ver_arn = partial(make_ver_arn, "layer")


def split_arn(arn_type: Callable[[str, str, str, str], str], arn: str) -> Any:
    for partition in PARTITION_NAMES:
        arn = arn.replace(f"arn:{partition}:lambda:", "")

    region, account, _, name, version = arn.split(":")

    return arn_type(region, account, name, version)


split_function_arn = partial(split_arn, ARN)
split_layer_arn = partial(split_arn, LAYER_ARN)


def get_backend(account_id: str, region: str) -> "LambdaBackend":
    from moto.core.models import default_user_config

    if default_user_config.get("lambda", {}).get("use_docker", True) is False:
        from moto.awslambda_simple.models import lambda_simple_backends

        return lambda_simple_backends[account_id][region]
    else:
        from moto.awslambda.models import lambda_backends

        return lambda_backends[account_id][region]
