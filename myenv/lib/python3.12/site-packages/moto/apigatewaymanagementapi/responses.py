"""Handles incoming apigatewaymanagementapi requests, invokes methods, returns responses."""

import json
from typing import Any

from moto.core.responses import BaseResponse

from .models import ApiGatewayManagementApiBackend, apigatewaymanagementapi_backends


class ApiGatewayManagementApiResponse(BaseResponse):
    """Handler for ApiGatewayManagementApi requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="apigatewaymanagementapi")

    def setup_class(
        self, request: Any, full_url: str, headers: Any, use_raw_body: bool = False
    ) -> None:
        super().setup_class(request, full_url, headers, use_raw_body=True)

    @property
    def apigatewaymanagementapi_backend(self) -> ApiGatewayManagementApiBackend:
        """Return backend instance specific for this region."""
        return apigatewaymanagementapi_backends[self.current_account][self.region]

    def delete_connection(self) -> str:
        connection_id = self.path.split("/@connections/")[-1]
        self.apigatewaymanagementapi_backend.delete_connection(
            connection_id=connection_id
        )
        return "{}"

    def get_connection(self) -> str:
        connection_id = self.path.split("/@connections/")[-1]
        connection = self.apigatewaymanagementapi_backend.get_connection(
            connection_id=connection_id
        )
        return json.dumps(connection.to_dict())

    def post_to_connection(self) -> str:
        connection_id = self.path.split("/@connections/")[-1]
        data = self.body
        self.apigatewaymanagementapi_backend.post_to_connection(
            data=data,
            connection_id=connection_id,
        )
        return "{}"
