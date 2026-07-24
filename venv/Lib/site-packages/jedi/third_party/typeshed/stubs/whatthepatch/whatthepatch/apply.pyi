from collections.abc import Iterable

from . import patch as patch
from .exceptions import HunkApplyException as HunkApplyException, SubprocessException as SubprocessException
from .snippets import remove as remove, which as which

def apply_patch(diffs: patch.diffobj | Iterable[patch.diffobj]) -> None: ...
def apply_diff(diff: patch.diffobj, text: str | Iterable[str], reverse: bool = False, use_patch: bool = False) -> list[str]: ...
