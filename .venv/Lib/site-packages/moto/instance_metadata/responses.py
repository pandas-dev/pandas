import datetime
import json
from typing import Any
from urllib.parse import urlparse

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse
from moto.core.utils import utcnow


class InstanceMetadataResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name=None)

    def backends(self) -> None:
        pass

    @staticmethod
    def metadata_response(  # type: ignore
        request: Any,  # pylint: disable=unused-argument
        full_url: str,
        headers: Any,
    ) -> TYPE_RESPONSE:
        """
        Mock response for localhost metadata

        http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/AESDG-chapter-instancedata.html
        """

        parsed_url = urlparse(full_url)
        tomorrow = utcnow() + datetime.timedelta(days=1)
        credentials = dict(
            AccessKeyId="test-key",
            SecretAccessKey="test-secret-key",
            Token="test-session-token",
            Expiration=tomorrow.strftime("%Y-%m-%dT%H:%M:%SZ"),
        )

        path = parsed_url.path

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
        else:
            raise NotImplementedError(
                f"The {path} metadata path has not been implemented"
            )
        try:
            from werkzeug.datastructures.headers import EnvironHeaders

            if isinstance(headers, EnvironHeaders):
                # We should be returning a generic dict here, not werkzeug-specific classes
                headers = {key: value for key, value in headers}
        except ImportError:
            pass
        return 200, headers, result
