from hvac.api.system_backend.system_backend_mixin import SystemBackendMixin

class Init(SystemBackendMixin):
    def read_init_status(self): ...
    def is_initialized(self): ...
    def initialize(
        self,
        secret_shares=None,
        secret_threshold=None,
        pgp_keys=None,
        root_token_pgp_key=None,
        stored_shares=None,
        recovery_shares=None,
        recovery_threshold=None,
        recovery_pgp_keys=None,
    ): ...
