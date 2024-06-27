import boto3
from botocore.client import BaseClient
from botocore.config import Config

from moto.core.models import botocore_stubber
from moto.settings import moto_server_port
from moto.stepfunctions.parser.asl.component.common.timeouts.timeout import (
    TimeoutSeconds,
)


def boto_client_for(region: str, account: str, service: str) -> BaseClient:
    intercepting_boto_calls = botocore_stubber.enabled
    kwargs = {}
    if not intercepting_boto_calls:
        kwargs["endpoint_url"] = f"http://localhost:{moto_server_port()}"
    return boto3.client(
        aws_access_key_id=account,
        region_name=region,
        service_name=service,
        aws_secret_access_key="sk",
        config=Config(
            parameter_validation=False,
            retries={"max_attempts": 0, "total_max_attempts": 1},
            connect_timeout=TimeoutSeconds.DEFAULT_TIMEOUT_SECONDS,
            read_timeout=TimeoutSeconds.DEFAULT_TIMEOUT_SECONDS,
        ),
        **kwargs,
    )
