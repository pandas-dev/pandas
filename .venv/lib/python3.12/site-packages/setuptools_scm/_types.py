from __future__ import annotations

import os

from typing import TYPE_CHECKING
from typing import Callable
from typing import List
from typing import Protocol
from typing import Sequence
from typing import Tuple
from typing import Union

from setuptools import Distribution

if TYPE_CHECKING:
    import sys

    if sys.version_info >= (3, 10):
        from typing import TypeAlias
    else:
        from typing_extensions import TypeAlias

    from . import version
    from ._integration.pyproject_reading import PyProjectData
    from ._integration.toml import InvalidTomlError

PathT: TypeAlias = Union["os.PathLike[str]", str]

CMD_TYPE: TypeAlias = Union[Sequence[PathT], str]

VERSION_SCHEME: TypeAlias = Union[str, Callable[["version.ScmVersion"], str]]
VERSION_SCHEMES: TypeAlias = Union[List[str], Tuple[str, ...], VERSION_SCHEME]
SCMVERSION: TypeAlias = "version.ScmVersion"

# Git pre-parse function types
GIT_PRE_PARSE: TypeAlias = Union[str, None]

# Testing injection types for configuration reading
GivenPyProjectResult: TypeAlias = Union[
    "PyProjectData", "InvalidTomlError", FileNotFoundError, None
]


class VersionInferenceApplicable(Protocol):
    """A result object from version inference decision that can be applied to a dist."""

    def apply(self, dist: Distribution) -> None:  # pragma: no cover - structural type
        ...


class GetVersionInferenceConfig(Protocol):
    """Callable protocol for the decision function used by integration points."""

    def __call__(
        self,
        dist_name: str | None,
        current_version: str | None,
        pyproject_data: PyProjectData,
        overrides: dict[str, object] | None = None,
    ) -> VersionInferenceApplicable:  # pragma: no cover - structural type
        ...
