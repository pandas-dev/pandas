from hvac.api.vault_api_base import VaultApiBase

DEFAULT_MOUNT_POINT: str

class Token(VaultApiBase):
    def create(
        self,
        id=None,
        role_name=None,
        policies=None,
        meta=None,
        no_parent: bool = False,
        no_default_policy: bool = False,
        renewable: bool = True,
        ttl=None,
        type=None,
        explicit_max_ttl=None,
        display_name: str = "token",
        num_uses: int = 0,
        period=None,
        entity_alias=None,
        wrap_ttl=None,
        mount_point="token",
    ): ...
    def create_orphan(
        self,
        id=None,
        role_name=None,
        policies=None,
        meta=None,
        no_default_policy: bool = False,
        renewable: bool = True,
        ttl=None,
        type=None,
        explicit_max_ttl=None,
        display_name: str = "token",
        num_uses: int = 0,
        period=None,
        entity_alias=None,
        wrap_ttl=None,
        mount_point="token",
    ): ...
    def list_accessors(self, mount_point="token"): ...
    def lookup(self, token, mount_point="token"): ...
    def lookup_self(self, mount_point="token"): ...
    def lookup_accessor(self, accessor, mount_point="token"): ...
    def renew(self, token, increment=None, wrap_ttl=None, mount_point="token"): ...
    def renew_self(self, increment=None, wrap_ttl=None, mount_point="token"): ...
    def renew_accessor(self, accessor, increment=None, wrap_ttl=None, mount_point="token"): ...
    def revoke(self, token, mount_point="token"): ...
    def revoke_self(self, mount_point="token"): ...
    def revoke_accessor(self, accessor, mount_point="token"): ...
    def revoke_and_orphan_children(self, token, mount_point="token"): ...
    def read_role(self, role_name, mount_point="token"): ...
    def list_roles(self, mount_point="token"): ...
    def create_or_update_role(
        self,
        role_name,
        allowed_policies=None,
        disallowed_policies=None,
        orphan: bool = False,
        renewable: bool = True,
        path_suffix=None,
        allowed_entity_aliases=None,
        mount_point="token",
        token_period=None,
        token_explicit_max_ttl=None,
    ): ...
    def delete_role(self, role_name, mount_point="token"): ...
    def tidy(self, mount_point="token"): ...
