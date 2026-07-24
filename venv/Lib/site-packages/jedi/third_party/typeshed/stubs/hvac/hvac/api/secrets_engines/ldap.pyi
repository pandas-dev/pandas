from typing import Final

from hvac.api.vault_api_base import VaultApiBase

DEFAULT_MOUNT_POINT: Final = "ldap"

class Ldap(VaultApiBase):
    def configure(
        self,
        binddn: str | None = None,
        bindpass: str | None = None,
        url: str | None = None,
        password_policy: str | None = None,
        schema: str | None = None,
        userdn: str | None = None,
        userattr: str | None = None,
        upndomain: str | None = None,
        connection_timeout: int | str | None = None,
        request_timeout: int | str | None = None,
        starttls: bool | None = None,
        insecure_tls: bool | None = None,
        certificate: str | None = None,
        client_tls_cert: str | None = None,
        client_tls_key: str | None = None,
        mount_point: str = "ldap",
    ): ...
    def read_config(self, mount_point: str = "ldap"): ...
    def rotate_root(self, mount_point: str = "ldap"): ...
    def create_or_update_static_role(
        self,
        name: str,
        username: str | None = None,
        dn: str | None = None,
        rotation_period: str | None = None,
        mount_point: str = "ldap",
    ): ...
    def read_static_role(self, name: str, mount_point: str = "ldap"): ...
    def list_static_roles(self, mount_point: str = "ldap"): ...
    def delete_static_role(self, name: str, mount_point: str = "ldap"): ...
    def generate_static_credentials(self, name: str, mount_point: str = "ldap"): ...
    def rotate_static_credentials(self, name: str, mount_point: str = "ldap"): ...
