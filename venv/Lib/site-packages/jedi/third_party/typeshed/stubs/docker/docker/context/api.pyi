from collections.abc import Sequence

from docker.context.context import Context
from docker.tls import TLSConfig

class ContextAPI:
    DEFAULT_CONTEXT: Context
    @classmethod
    def create_context(
        cls,
        name: str,
        orchestrator: str | None = None,
        host: str | None = None,
        tls_cfg: TLSConfig | None = None,
        default_namespace: str | None = None,
        skip_tls_verify: bool = False,
    ) -> Context: ...
    @classmethod
    def get_context(cls, name: str | None = None) -> Context: ...
    @classmethod
    def contexts(cls) -> Sequence[Context]: ...
    @classmethod
    def get_current_context(cls) -> Context: ...
    @classmethod
    def set_current_context(cls, name: str = "default") -> None: ...
    @classmethod
    def remove_context(cls, name: str) -> None: ...
    @classmethod
    def inspect_context(cls, name: str = "default") -> Context: ...
