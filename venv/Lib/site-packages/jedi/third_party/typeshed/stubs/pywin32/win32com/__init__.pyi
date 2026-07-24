from collections.abc import MutableSequence

from . import gen_py as gen_py

__gen_path__: str
__build_path__: str | None

def SetupEnvironment() -> None: ...
def __PackageSupportBuildPath__(package_path: MutableSequence[str]) -> None: ...
