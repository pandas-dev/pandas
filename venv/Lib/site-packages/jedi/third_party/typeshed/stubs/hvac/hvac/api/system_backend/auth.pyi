from hvac.api.system_backend.system_backend_mixin import SystemBackendMixin

class Auth(SystemBackendMixin):
    def list_auth_methods(self): ...
    def enable_auth_method(
        self, method_type, description=None, config=None, plugin_name=None, local: bool = False, path=None, **kwargs
    ): ...
    def disable_auth_method(self, path): ...
    def read_auth_method_tuning(self, path): ...
    def tune_auth_method(
        self,
        path,
        default_lease_ttl=None,
        max_lease_ttl=None,
        description=None,
        audit_non_hmac_request_keys=None,
        audit_non_hmac_response_keys=None,
        listing_visibility=None,
        passthrough_request_headers=None,
        **kwargs,
    ): ...
