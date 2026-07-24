from hvac.api.vault_api_base import VaultApiBase

DEFAULT_MOUNT_POINT: str

class ActiveDirectory(VaultApiBase):
    def configure(
        self,
        binddn=None,
        bindpass=None,
        url=None,
        userdn=None,
        upndomain=None,
        ttl=None,
        max_ttl=None,
        mount_point="ad",
        *args,
        **kwargs,
    ): ...
    def read_config(self, mount_point="ad"): ...
    def create_or_update_role(self, name, service_account_name=None, ttl=None, mount_point="ad"): ...
    def read_role(self, name, mount_point="ad"): ...
    def list_roles(self, mount_point="ad"): ...
    def delete_role(self, name, mount_point="ad"): ...
    def generate_credentials(self, name, mount_point="ad"): ...
