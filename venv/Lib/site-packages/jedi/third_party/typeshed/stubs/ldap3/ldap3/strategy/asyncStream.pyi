from _typeshed import Incomplete

from ..strategy.asynchronous import AsyncStrategy

class AsyncStreamStrategy(AsyncStrategy):
    can_stream: bool
    line_separator: Incomplete
    all_base64: bool
    stream: Incomplete
    order: Incomplete
    persistent_search_message_id: Incomplete
    streaming: bool
    callback: Incomplete
    events: Incomplete
    def __init__(self, ldap_connection) -> None: ...
    def accumulate_stream(self, message_id, change) -> None: ...
    def get_stream(self): ...
    def set_stream(self, value) -> None: ...
