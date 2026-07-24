from typing import Any

from gevent.hub import Hub
from gevent.resolver import AbstractResolver

class Resolver(AbstractResolver):
    def __init__(self, hub: Hub | None = ...) -> None: ...
    @property
    def resolver(self) -> Any: ...  # this is a custom dnspython Resolver

__all__ = ["Resolver"]
