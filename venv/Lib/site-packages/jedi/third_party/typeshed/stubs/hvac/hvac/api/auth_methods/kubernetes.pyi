from hvac.api.vault_api_base import VaultApiBase

DEFAULT_MOUNT_POINT: str

class Kubernetes(VaultApiBase):
    def configure(
        self,
        kubernetes_host,
        kubernetes_ca_cert=None,
        token_reviewer_jwt=None,
        pem_keys=None,
        issuer=None,
        mount_point="kubernetes",
        disable_local_ca_jwt: bool = False,
    ): ...
    def read_config(self, mount_point="kubernetes"): ...
    def create_role(
        self,
        name,
        bound_service_account_names,
        bound_service_account_namespaces,
        ttl=None,
        max_ttl=None,
        period=None,
        policies=None,
        token_type: str = "",
        mount_point="kubernetes",
        alias_name_source=None,
        audience: str | None = None,
    ): ...
    def read_role(self, name, mount_point="kubernetes"): ...
    def list_roles(self, mount_point="kubernetes"): ...
    def delete_role(self, name, mount_point="kubernetes"): ...
    def login(self, role, jwt, use_token: bool = True, mount_point="kubernetes"): ...
