from hvac.api.system_backend.system_backend_mixin import SystemBackendMixin

class Capabilities(SystemBackendMixin):
    def get_capabilities(self, paths, token=None, accessor=None): ...
