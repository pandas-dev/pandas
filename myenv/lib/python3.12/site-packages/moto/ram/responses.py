import json
from typing import Any, Dict

from moto.core.responses import BaseResponse

from .models import ResourceAccessManagerBackend, ram_backends


class ResourceAccessManagerResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="ram")

    @property
    def ram_backend(self) -> ResourceAccessManagerBackend:
        return ram_backends[self.current_account][self.region]

    @property
    def request_params(self) -> Dict[str, Any]:  # type: ignore[misc]
        try:
            return json.loads(self.body)
        except ValueError:
            return {}

    def create_resource_share(self) -> str:
        return json.dumps(self.ram_backend.create_resource_share(**self.request_params))

    def get_resource_shares(self) -> str:
        return json.dumps(self.ram_backend.get_resource_shares(**self.request_params))

    def update_resource_share(self) -> str:
        return json.dumps(self.ram_backend.update_resource_share(**self.request_params))

    def delete_resource_share(self) -> str:
        arn = self._get_param("resourceShareArn")
        return json.dumps(self.ram_backend.delete_resource_share(arn))

    def enable_sharing_with_aws_organization(self) -> str:
        return json.dumps(self.ram_backend.enable_sharing_with_aws_organization())
