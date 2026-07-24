from hvac.api.system_backend.system_backend_mixin import SystemBackendMixin

class Health(SystemBackendMixin):
    def read_health_status(
        self,
        standby_ok=None,
        active_code=None,
        standby_code=None,
        dr_secondary_code=None,
        performance_standby_code=None,
        sealed_code=None,
        uninit_code=None,
        method: str = "HEAD",
    ): ...
