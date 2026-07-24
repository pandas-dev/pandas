from _typeshed import Incomplete

class ConnectionUsage:
    open_sockets: int
    closed_sockets: int
    wrapped_sockets: int
    bytes_transmitted: int
    bytes_received: int
    messages_transmitted: int
    messages_received: int
    operations: int
    abandon_operations: int
    add_operations: int
    bind_operations: int
    compare_operations: int
    delete_operations: int
    extended_operations: int
    modify_operations: int
    modify_dn_operations: int
    search_operations: int
    unbind_operations: int
    referrals_received: int
    referrals_followed: int
    referrals_connections: int
    restartable_failures: int
    restartable_successes: int
    servers_from_pool: int
    def reset(self) -> None: ...
    initial_connection_start_time: Incomplete
    open_socket_start_time: Incomplete
    connection_stop_time: Incomplete
    last_transmitted_time: Incomplete
    last_received_time: Incomplete
    def __init__(self) -> None: ...
    def __iadd__(self, other): ...
    def update_transmitted_message(self, message, length) -> None: ...
    def update_received_message(self, length) -> None: ...
    def start(self, reset: bool = True) -> None: ...
    def stop(self) -> None: ...
    @property
    def elapsed_time(self): ...
