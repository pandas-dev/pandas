from hvac.api.vault_api_base import VaultApiBase

DEFAULT_MOUNT_POINT: str

class Ldap(VaultApiBase):
    def configure(
        self,
        userdn=None,
        groupdn=None,
        url=None,
        case_sensitive_names=None,
        starttls=None,
        tls_min_version=None,
        tls_max_version=None,
        insecure_tls=None,
        certificate=None,
        binddn=None,
        bindpass=None,
        userattr=None,
        discoverdn=None,
        deny_null_bind: bool = True,
        upndomain=None,
        groupfilter=None,
        groupattr=None,
        use_token_groups=None,
        token_ttl=None,
        token_max_ttl=None,
        mount_point="ldap",
        *,
        anonymous_group_search=None,
        client_tls_cert=None,
        client_tls_key=None,
        connection_timeout=None,
        dereference_aliases=None,
        max_page_size=None,
        request_timeout=None,
        token_bound_cidrs=None,
        token_explicit_max_ttl=None,
        token_no_default_policy=None,
        token_num_uses=None,
        token_period=None,
        token_policies=None,
        token_type=None,
        userfilter=None,
        username_as_alias=None,
    ): ...
    def read_configuration(self, mount_point="ldap"): ...
    def create_or_update_group(self, name, policies=None, mount_point="ldap"): ...
    def list_groups(self, mount_point="ldap"): ...
    def read_group(self, name, mount_point="ldap"): ...
    def delete_group(self, name, mount_point="ldap"): ...
    def create_or_update_user(self, username, policies=None, groups=None, mount_point="ldap"): ...
    def list_users(self, mount_point="ldap"): ...
    def read_user(self, username, mount_point="ldap"): ...
    def delete_user(self, username, mount_point="ldap"): ...
    def login(self, username, password, use_token: bool = True, mount_point="ldap"): ...
