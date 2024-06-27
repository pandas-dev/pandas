from typing import Tuple, Union

import requests

from ..models import Integration
from . import IntegrationParser


class TypeHttpParser(IntegrationParser):
    """
    Parse invocations to a APIGateway resource with integration type HTTP
    """

    def invoke(
        self, request: requests.PreparedRequest, integration: Integration
    ) -> Tuple[int, Union[str, bytes]]:
        uri = integration.uri
        requests_func = getattr(requests, integration.http_method.lower())  # type: ignore[union-attr]
        response = requests_func(uri)
        return response.status_code, response.text
