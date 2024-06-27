import re
from io import BytesIO
from typing import Any, Optional, Union

from botocore.awsrequest import AWSResponse

import moto.backend_index as backend_index
from moto.core.base_backend import BackendDict
from moto.core.common_types import TYPE_RESPONSE
from moto.core.config import passthrough_service, passthrough_url
from moto.core.utils import get_equivalent_url_in_aws_domain


class MockRawResponse(BytesIO):
    def __init__(self, response_input: Union[str, bytes]):
        if isinstance(response_input, str):
            response_input = response_input.encode("utf-8")
        super().__init__(response_input)

    def stream(self, **kwargs: Any) -> Any:  # pylint: disable=unused-argument
        contents = self.read()
        while contents:
            yield contents
            contents = self.read()


class BotocoreStubber:
    def __init__(self) -> None:
        self.enabled = False

    def __call__(
        self, event_name: str, request: Any, **kwargs: Any
    ) -> Optional[AWSResponse]:
        if not self.enabled:
            return None

        response = self.process_request(request)
        if response is not None:
            status, headers, body = response
            return AWSResponse(request.url, status, headers, MockRawResponse(body))  # type: ignore[arg-type]
        else:
            return response

    def process_request(self, request: Any) -> Optional[TYPE_RESPONSE]:
        # Handle non-standard AWS endpoint hostnames from ISO regions or custom
        # S3 endpoints.
        parsed_url, _ = get_equivalent_url_in_aws_domain(request.url)
        # Remove the querystring from the URL, as we'll never match on that
        clean_url = f"{parsed_url.scheme}://{parsed_url.netloc}{parsed_url.path}"

        if passthrough_url(clean_url):
            return None

        for service, pattern in backend_index.backend_url_patterns:
            if pattern.match(clean_url):
                if passthrough_service(service):
                    return None

                import moto.backends as backends
                from moto.core import DEFAULT_ACCOUNT_ID
                from moto.core.exceptions import HTTPException

                # TODO: cache this part - we only need backend.urls
                backend_dict = backends.get_backend(service)  # type: ignore[call-overload]

                if isinstance(backend_dict, BackendDict):
                    if "us-east-1" in backend_dict[DEFAULT_ACCOUNT_ID]:
                        backend = backend_dict[DEFAULT_ACCOUNT_ID]["us-east-1"]
                    else:
                        backend = backend_dict[DEFAULT_ACCOUNT_ID]["aws"]
                else:
                    backend = backend_dict["global"]

                for header, value in request.headers.items():
                    if isinstance(value, bytes):
                        request.headers[header] = value.decode("utf-8")

                for url, method_to_execute in backend.urls.items():
                    if re.compile(url).match(clean_url):
                        from moto.moto_api import recorder

                        try:
                            recorder._record_request(request)
                            status, headers, body = method_to_execute(
                                request, request.url, request.headers
                            )
                        except HTTPException as e:
                            status = e.code
                            headers = e.get_headers()
                            body = e.get_body()

                        return status, headers, body

        if re.compile(r"https?://.+\.amazonaws.com/.*").match(clean_url):
            return 404, {}, "Not yet implemented"

        return None
