from .restartable import RestartableStrategy

class SafeRestartableStrategy(RestartableStrategy):
    thread_safe: bool
    def __init__(self, ldap_connection) -> None: ...
