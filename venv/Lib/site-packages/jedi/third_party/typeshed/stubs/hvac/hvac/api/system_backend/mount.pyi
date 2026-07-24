from hvac.api.system_backend.system_backend_mixin import SystemBackendMixin

class Mount(SystemBackendMixin):
    def list_mounted_secrets_engines(self): ...
    def retrieve_mount_option(self, mount_point, option_name, default_value=None): ...
    def enable_secrets_engine(
        self,
        backend_type,
        path=None,
        description=None,
        config=None,
        plugin_name=None,
        options=None,
        local: bool = False,
        seal_wrap: bool = False,
        **kwargs,
    ): ...
    def disable_secrets_engine(self, path): ...
    def read_mount_configuration(self, path): ...
    def tune_mount_configuration(
        self,
        path,
        default_lease_ttl=None,
        max_lease_ttl=None,
        description=None,
        audit_non_hmac_request_keys=None,
        audit_non_hmac_response_keys=None,
        listing_visibility=None,
        passthrough_request_headers=None,
        options=None,
        force_no_cache=None,
        **kwargs,
    ): ...
    def move_backend(self, from_path, to_path): ...
