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
        resource_owner = self._get_param("resourceOwner")
        return json.dumps(self.ram_backend.get_resource_shares(resource_owner))

    def update_resource_share(self) -> str:
        resource_share_arn = self._get_param("resourceShareArn")
        allow_external_principals = self._get_param("allowExternalPrincipals", True)
        name = self._get_param("name")
        return json.dumps(
            self.ram_backend.update_resource_share(
                resource_share_arn, allow_external_principals, name
            )
        )

    def delete_resource_share(self) -> str:
        resource_share_arn = self._get_param("resourceShareArn")
        return json.dumps(self.ram_backend.delete_resource_share(resource_share_arn))

    def enable_sharing_with_aws_organization(self) -> str:
        return json.dumps(self.ram_backend.enable_sharing_with_aws_organization())

    def get_resource_share_associations(self) -> str:
        association_type = self._get_param("associationType")
        association_status = self._get_param("associationStatus")
        resource_share_arns = self._get_param("resourceShareArns", [])
        resource_arn = self._get_param("resourceArn")
        principal = self._get_param("principal")
        return json.dumps(
            self.ram_backend.get_resource_share_associations(
                association_type,
                association_status,
                resource_share_arns,
                resource_arn,
                principal,
            )
        )

    def list_resource_types(self) -> str:
        resource_region_scope = self._get_param("resourceRegionScope", "ALL")
        return json.dumps(self.ram_backend.list_resource_types(resource_region_scope))

    def list_permissions(self) -> str:
        resource_type = self._get_param("resourceType")
        permission_type = self._get_param("permissionType", "ALL")
        return json.dumps(
            self.ram_backend.list_permissions(resource_type, permission_type)
        )
