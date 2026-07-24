from .sync import SyncStrategy

class SafeSyncStrategy(SyncStrategy):
    thread_safe: bool
    def __init__(self, ldap_connection) -> None: ...
