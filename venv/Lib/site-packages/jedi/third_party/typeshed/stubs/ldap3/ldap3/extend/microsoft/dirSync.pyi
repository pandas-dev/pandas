from _typeshed import Incomplete

class DirSync:
    connection: Incomplete
    base: Incomplete
    filter: Incomplete
    attributes: Incomplete
    cookie: Incomplete
    object_security: Incomplete
    ancestors_first: Incomplete
    public_data_only: Incomplete
    incremental_values: Incomplete
    max_length: Incomplete
    hex_guid: Incomplete
    more_results: bool
    def __init__(
        self,
        connection,
        sync_base,
        sync_filter,
        attributes,
        cookie,
        object_security,
        ancestors_first,
        public_data_only,
        incremental_values,
        max_length,
        hex_guid,
    ) -> None: ...
    def loop(self): ...
