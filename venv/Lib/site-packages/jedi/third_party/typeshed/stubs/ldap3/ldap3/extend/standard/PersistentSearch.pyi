from _typeshed import Incomplete

class PersistentSearch:
    connection: Incomplete
    changes_only: Incomplete
    notifications: Incomplete
    message_id: Incomplete
    base: Incomplete
    filter: Incomplete
    scope: Incomplete
    dereference_aliases: Incomplete
    attributes: Incomplete
    size_limit: Incomplete
    time_limit: Incomplete
    controls: Incomplete
    def __init__(
        self,
        connection,
        search_base,
        search_filter,
        search_scope,
        dereference_aliases,
        attributes,
        size_limit,
        time_limit,
        controls,
        changes_only,
        events_type,
        notifications,
        streaming,
        callback,
    ) -> None: ...
    def start(self) -> None: ...
    def stop(self, unbind: bool = True) -> None: ...
    def next(self, block: bool = False, timeout=None): ...
    def funnel(self, block: bool = False, timeout=None) -> None: ...
