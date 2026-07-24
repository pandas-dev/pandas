from _typeshed import Unused
from collections.abc import Mapping, Sequence

def launch(
    *,
    program: str | None = None,
    cwd: str | None = None,
    args: Sequence[str] = (),
    env: Mapping[str, str] | None = None,
    stopAtBeginningOfMainSubprogram: bool = False,
    **extra: Unused,
): ...
def attach(*, program: str | None = None, pid: int | None = None, target: str | None = None, **args: Unused) -> None: ...
def config_done(**args: Unused) -> None: ...
