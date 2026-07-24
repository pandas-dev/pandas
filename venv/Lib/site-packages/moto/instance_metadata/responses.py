import datetime
import json
from typing import Any
from urllib.parse import urlparse

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse
from moto.core.utils import utcnow

INSTANCE_METADATA = {
    "ami-id": "ami-moto-virtual",
    "ami-launch-index": "0",
    "ami-manifest-path": "(unknown)",
    "hostname": "moto-instance.ec2.internal",
    "instance-action": "none",
    "instance-id": "i-fake",
    "instance-life-cycle": "on-demand",
    "instance-type": "moto.virtual",
    "local-hostname": "moto-instance.ec2.internal",
    "local-ipv4": "127.0.0.1",
    "mac": "00:de:ad:be:ef:00",
    "placement/availability-zone": "us-east-1a",
    "placement/availability-zone-id": "use1-az1",
    "placement/region": "us-east-1",
    "profile": "default-hvm",
    "public-hostname": "moto-instance.compute-1.amazonaws.com",
    "public-ipv4": "198.51.100.1",
    "reservation-id": "r-moto-virtual",
    "security-groups": "moto-default",
}


class InstanceMetadataResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name=None)

    def backends(self) -> None:
        pass

    @staticmethod
    def metadata_response(  # type: ignore
        request: Any,
        full_url: str,
        headers: Any,
    ) -> TYPE_RESPONSE:
        """
        Mock response for localhost metadata

        http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/AESDG-chapter-instancedata.html
        """

        parsed_url = urlparse(full_url)
        tomorrow = utcnow() + datetime.timedelta(days=1)
        credentials = {
            "AccessKeyId": "test-key",
            "SecretAccessKey": "test-secret-key",
            "Token": "test-session-token",
            "Expiration": tomorrow.strftime("%Y-%m-%dT%H:%M:%SZ"),
        }

        path = parsed_url.path

        # IMDSv2 token endpoint
        if path == "/latest/api/token":
            result = "fake-imdsv2-token"
            return 200, {}, result

        meta_data_prefix = "/latest/meta-data/"
        # Strip prefix if it is there
        if path.startswith(meta_data_prefix):
            path = path[len(meta_data_prefix) :]

        if path == "":
            result = "iam"
        elif path == "iam":
            result = json.dumps({"security-credentials": {"default-role": credentials}})
        elif path == "iam/security-credentials/":
            result = "default-role"
        elif path == "iam/security-credentials/default-role":
            result = json.dumps(credentials)
        elif path in INSTANCE_METADATA:
            result = INSTANCE_METADATA[path]
        else:
            raise NotImplementedError(
                f"The {path} metadata path has not been implemented"
            )
        try:
            from werkzeug.datastructures.headers import EnvironHeaders

            if isinstance(headers, EnvironHeaders):
                # We should be returning a generic dict here, not werkzeug-specific classes
                headers = dict(headers)
        except ImportError:
            pass
        return 200, headers, result
