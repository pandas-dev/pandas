"""configuration"""

from __future__ import annotations

import dataclasses
import os
import re
import warnings

from pathlib import Path
from typing import TYPE_CHECKING
from typing import Any
from typing import Pattern
from typing import Protocol

if TYPE_CHECKING:
    from . import git

from . import _log
from . import _types as _t
from ._integration.pyproject_reading import PyProjectData
from ._integration.pyproject_reading import (
    get_args_for_pyproject as _get_args_for_pyproject,
)
from ._integration.pyproject_reading import read_pyproject as _read_pyproject
from ._overrides import read_toml_overrides
from ._version_cls import Version as _Version
from ._version_cls import _validate_version_cls
from ._version_cls import _VersionT

log = _log.log.getChild("config")


def _is_called_from_dataclasses() -> bool:
    """Check if the current call is from the dataclasses module."""
    import inspect

    frame = inspect.currentframe()
    try:
        # Walk up to 7 frames to check for dataclasses calls
        current_frame = frame
        assert current_frame is not None
        for _ in range(7):
            current_frame = current_frame.f_back
            if current_frame is None:
                break
            if "dataclasses.py" in current_frame.f_code.co_filename:
                return True
        return False
    finally:
        del frame


class _GitDescribeCommandDescriptor:
    """Data descriptor for deprecated git_describe_command field."""

    def __get__(
        self, obj: Configuration | None, objtype: type[Configuration] | None = None
    ) -> _t.CMD_TYPE | None:
        if obj is None:
            return self  # type: ignore[return-value]

        # Only warn if not being called by dataclasses.replace or similar introspection
        is_from_dataclasses = _is_called_from_dataclasses()
        if not is_from_dataclasses:
            warnings.warn(
                "Configuration field 'git_describe_command' is deprecated. "
                "Use 'scm.git.describe_command' instead.",
                DeprecationWarning,
                stacklevel=2,
            )
        return obj.scm.git.describe_command

    def __set__(self, obj: Configuration, value: _t.CMD_TYPE | None) -> None:
        warnings.warn(
            "Configuration field 'git_describe_command' is deprecated. "
            "Use 'scm.git.describe_command' instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        obj.scm.git.describe_command = value


DEFAULT_TAG_REGEX = re.compile(
    r"^(?:[\w-]+-)?(?P<version>[vV]?\d+(?:\.\d+){0,2}[^\+]*)(?:\+.*)?$"
)
"""default tag regex that tries to match PEP440 style versions
with prefix consisting of dashed words"""

DEFAULT_VERSION_SCHEME = "guess-next-dev"
DEFAULT_LOCAL_SCHEME = "node-and-date"


def _check_tag_regex(value: str | Pattern[str] | None) -> Pattern[str]:
    if not value:
        regex = DEFAULT_TAG_REGEX
    else:
        regex = re.compile(value)

    group_names = regex.groupindex.keys()
    if regex.groups == 0 or (regex.groups > 1 and "version" not in group_names):
        raise ValueError(
            f"Expected tag_regex '{regex.pattern}' to contain a single match group or"
            " a group named 'version' to identify the version part of any tag."
        )

    return regex


def _get_default_git_pre_parse() -> git.GitPreParse:
    """Get the default git pre_parse enum value"""
    from . import git

    return git.GitPreParse.WARN_ON_SHALLOW


class ParseFunction(Protocol):
    def __call__(
        self, root: _t.PathT, *, config: Configuration
    ) -> _t.SCMVERSION | None: ...


def _check_absolute_root(root: _t.PathT, relative_to: _t.PathT | None) -> str:
    log.debug("check absolute root=%s relative_to=%s", root, relative_to)
    if relative_to:
        if (
            os.path.isabs(root)
            and os.path.isabs(relative_to)
            and not os.path.commonpath([root, relative_to]) == root
        ):
            warnings.warn(
                f"absolute root path '{root}' overrides relative_to '{relative_to}'"
            )
        if os.path.isdir(relative_to):
            warnings.warn(
                "relative_to is expected to be a file,"
                f" its the directory {relative_to}\n"
                "assuming the parent directory was passed"
            )
            log.debug("dir %s", relative_to)
            root = os.path.join(relative_to, root)
        else:
            log.debug("file %s", relative_to)
            root = os.path.join(os.path.dirname(relative_to), root)
    return os.path.abspath(root)


@dataclasses.dataclass
class GitConfiguration:
    """Git-specific configuration options"""

    pre_parse: git.GitPreParse = dataclasses.field(
        default_factory=lambda: _get_default_git_pre_parse()
    )
    describe_command: _t.CMD_TYPE | None = None

    @classmethod
    def from_data(cls, data: dict[str, Any]) -> GitConfiguration:
        """Create GitConfiguration from configuration data, converting strings to enums"""
        git_data = data.copy()

        # Convert string pre_parse values to enum instances
        if "pre_parse" in git_data and isinstance(git_data["pre_parse"], str):
            from . import git

            try:
                git_data["pre_parse"] = git.GitPreParse(git_data["pre_parse"])
            except ValueError as e:
                valid_options = [option.value for option in git.GitPreParse]
                raise ValueError(
                    f"Invalid git pre_parse function '{git_data['pre_parse']}'. "
                    f"Valid options are: {', '.join(valid_options)}"
                ) from e

        return cls(**git_data)


@dataclasses.dataclass
class ScmConfiguration:
    """SCM-specific configuration options"""

    git: GitConfiguration = dataclasses.field(default_factory=GitConfiguration)

    @classmethod
    def from_data(cls, data: dict[str, Any]) -> ScmConfiguration:
        """Create ScmConfiguration from configuration data"""
        scm_data = data.copy()

        # Handle git-specific configuration
        git_data = scm_data.pop("git", {})
        git_config = GitConfiguration.from_data(git_data)

        return cls(git=git_config, **scm_data)


@dataclasses.dataclass
class Configuration:
    """Global configuration model"""

    relative_to: _t.PathT | None = None
    root: _t.PathT = "."
    version_scheme: _t.VERSION_SCHEME = DEFAULT_VERSION_SCHEME
    local_scheme: _t.VERSION_SCHEME = DEFAULT_LOCAL_SCHEME
    tag_regex: Pattern[str] = DEFAULT_TAG_REGEX
    parentdir_prefix_version: str | None = None
    fallback_version: str | None = None
    fallback_root: _t.PathT = "."
    write_to: _t.PathT | None = None
    write_to_template: str | None = None
    version_file: _t.PathT | None = None
    version_file_template: str | None = None
    parse: ParseFunction | None = None
    git_describe_command: dataclasses.InitVar[_t.CMD_TYPE | None] = (
        _GitDescribeCommandDescriptor()
    )

    dist_name: str | None = None
    version_cls: type[_VersionT] = _Version
    search_parent_directories: bool = False

    parent: _t.PathT | None = None

    # Nested SCM configurations
    scm: ScmConfiguration = dataclasses.field(
        default_factory=lambda: ScmConfiguration()
    )

    # Deprecated fields (handled in __post_init__)

    def __post_init__(self, git_describe_command: _t.CMD_TYPE | None) -> None:
        self.tag_regex = _check_tag_regex(self.tag_regex)

        # Handle deprecated git_describe_command
        # Check if it's a descriptor object (happens when no value is passed)
        if git_describe_command is not None and not isinstance(
            git_describe_command, _GitDescribeCommandDescriptor
        ):
            # Check if this is being called from dataclasses
            is_from_dataclasses = _is_called_from_dataclasses()

            same_value = (
                self.scm.git.describe_command is not None
                and self.scm.git.describe_command == git_describe_command
            )

            if is_from_dataclasses and same_value:
                # Ignore the passed value - it's from dataclasses.replace() with same value
                pass
            else:
                warnings.warn(
                    "Configuration field 'git_describe_command' is deprecated. "
                    "Use 'scm.git.describe_command' instead.",
                    DeprecationWarning,
                    stacklevel=2,
                )
                # Check for conflicts
                if self.scm.git.describe_command is not None:
                    raise ValueError(
                        "Cannot specify both 'git_describe_command' (deprecated) and "
                        "'scm.git.describe_command'. Please use only 'scm.git.describe_command'."
                    )
                self.scm.git.describe_command = git_describe_command

    @property
    def absolute_root(self) -> str:
        return _check_absolute_root(self.root, self.relative_to)

    @classmethod
    def from_file(
        cls,
        name: str | os.PathLike[str] = "pyproject.toml",
        dist_name: str | None = None,
        pyproject_data: PyProjectData | None = None,
        **kwargs: Any,
    ) -> Configuration:
        """
                Read Configuration from pyproject.toml (or similar).
                Raises exceptions when file is not found or toml is
                not installed or the file has invalid format.

        Parameters:
        - name: path to pyproject.toml
        - dist_name: name of the distribution
        - **kwargs: additional keyword arguments to pass to the Configuration constructor
        """

        if pyproject_data is None:
            pyproject_data = _read_pyproject(Path(name))
        args = _get_args_for_pyproject(pyproject_data, dist_name, kwargs)

        args.update(read_toml_overrides(args["dist_name"]))
        relative_to = args.pop("relative_to", name)
        return cls.from_data(relative_to=relative_to, data=args)

    @classmethod
    def from_data(
        cls, relative_to: str | os.PathLike[str], data: dict[str, Any]
    ) -> Configuration:
        """
        given configuration data
        create a config instance after validating tag regex/version class
        """
        version_cls = _validate_version_cls(
            data.pop("version_cls", None), data.pop("normalize", True)
        )

        # Handle nested SCM configuration
        scm_data = data.pop("scm", {})

        # Handle nested SCM configuration

        scm_config = ScmConfiguration.from_data(scm_data)
        return cls(
            relative_to=relative_to,
            version_cls=version_cls,
            scm=scm_config,
            **data,
        )
