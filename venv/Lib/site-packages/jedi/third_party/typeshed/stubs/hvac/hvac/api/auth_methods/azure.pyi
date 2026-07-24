import logging

from hvac.api.vault_api_base import VaultApiBase

DEFAULT_MOUNT_POINT: str

logger: logging.Logger

class Azure(VaultApiBase):
    def configure(self, tenant_id, resource, environment=None, client_id=None, client_secret=None, mount_point="azure"): ...
    def read_config(self, mount_point="azure"): ...
    def delete_config(self, mount_point="azure"): ...
    def create_role(
        self,
        name,
        policies=None,
        ttl=None,
        max_ttl=None,
        period=None,
        bound_service_principal_ids=None,
        bound_group_ids=None,
        bound_locations=None,
        bound_subscription_ids=None,
        bound_resource_groups=None,
        bound_scale_sets=None,
        num_uses=None,
        mount_point="azure",
    ): ...
    def read_role(self, name, mount_point="azure"): ...
    def list_roles(self, mount_point="azure"): ...
    def delete_role(self, name, mount_point="azure"): ...
    def login(
        self,
        role,
        jwt,
        subscription_id=None,
        resource_group_name=None,
        vm_name=None,
        vmss_name=None,
        use_token: bool = True,
        mount_point="azure",
    ): ...
