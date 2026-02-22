from __future__ import annotations

import dataclasses
import logging
import re
import warnings

from pathlib import Path
from typing import Any
from typing import NoReturn
from typing import Pattern

from . import _config
from . import _entrypoints
from . import _run_cmd
from . import _types as _t
from ._config import Configuration
from ._overrides import _read_pretended_version_for
from ._version_cls import _validate_version_cls
from .version import ScmVersion
from .version import format_version as _format_version

EMPTY_TAG_REGEX_DEPRECATION = DeprecationWarning(
    "empty regex for tag regex is invalid, using default"
)

_log = logging.getLogger(__name__)


def parse_scm_version(config: Configuration) -> ScmVersion | None:
    try:
        if config.parse is not None:
            parse_result = config.parse(config.absolute_root, config=config)
            if parse_result is not None and not isinstance(parse_result, ScmVersion):
                raise TypeError(
                    f"version parse result was {str!r}\n"
                    "please return a parsed version (ScmVersion)"
                )
            return parse_result
        else:
            return _entrypoints.version_from_entrypoint(
                config,
                entrypoint="setuptools_scm.parse_scm",
                root=config.absolute_root,
            )
    except _run_cmd.CommandNotFoundError as e:
        _log.exception("command %s not found while parsing the scm, using fallbacks", e)
        return None


def parse_fallback_version(config: Configuration) -> ScmVersion | None:
    return _entrypoints.version_from_entrypoint(
        config,
        entrypoint="setuptools_scm.parse_scm_fallback",
        root=config.fallback_root,
    )


def parse_version(config: Configuration) -> ScmVersion | None:
    # First try to get a version from the normal flow
    scm_version = (
        _read_pretended_version_for(config)
        or parse_scm_version(config)
        or parse_fallback_version(config)
    )

    # Apply any metadata overrides to the version we found
    from ._overrides import _apply_metadata_overrides

    return _apply_metadata_overrides(scm_version, config)


def write_version_files(
    config: Configuration, version: str, scm_version: ScmVersion
) -> None:
    if config.write_to is not None:
        from ._integration.dump_version import dump_version

        dump_version(
            root=config.root,
            version=version,
            scm_version=scm_version,
            write_to=config.write_to,
            template=config.write_to_template,
        )
    if config.version_file:
        from ._integration.dump_version import write_version_to_path

        version_file = Path(config.version_file)
        assert not version_file.is_absolute(), f"{version_file=}"
        # todo: use a better name than fallback root
        assert config.relative_to is not None
        target = Path(config.relative_to).parent.joinpath(version_file)
        write_version_to_path(
            target,
            template=config.version_file_template,
            version=version,
            scm_version=scm_version,
        )


def _get_version(
    config: Configuration, force_write_version_files: bool | None = None
) -> str | None:
    parsed_version = parse_version(config)
    if parsed_version is None:
        return None
    version_string = _format_version(parsed_version)
    if force_write_version_files is None:
        force_write_version_files = True
        warnings.warn(
            "force_write_version_files ought to be set,"
            " presuming the legacy True value",
            DeprecationWarning,
        )

    if force_write_version_files:
        write_version_files(config, version=version_string, scm_version=parsed_version)

    return version_string


def _find_scm_in_parents(config: Configuration) -> Path | None:
    """
    Search parent directories for SCM repositories when relative_to is not set.
    Uses the existing entrypoint system for SCM discovery.
    """
    if config.search_parent_directories:
        return None

    searching_config = dataclasses.replace(config, search_parent_directories=True)

    from .discover import iter_matching_entrypoints

    for _ep in iter_matching_entrypoints(
        config.absolute_root, "setuptools_scm.parse_scm", searching_config
    ):
        # xxx: iter_matching_entrypoints should return the parent directory, we do a hack atm
        assert searching_config.parent is not None
        return Path(searching_config.parent)

    return None


def _version_missing(config: Configuration) -> NoReturn:
    base_error = (
        f"setuptools-scm was unable to detect version for {config.absolute_root}.\n\n"
    )

    # If relative_to is not set, check for SCM repositories in parent directories
    scm_parent = None
    if config.relative_to is None:
        scm_parent = _find_scm_in_parents(config)

    if scm_parent is not None:
        # Found an SCM repository in a parent directory
        error_msg = (
            base_error
            + f"However, a repository was found in a parent directory: {scm_parent}\n\n"
            f"To fix this, you have a few options:\n\n"
            f"1. Use the 'relative_to' parameter to specify the file that setuptools-scm should use as reference:\n"
            f"   setuptools_scm.get_version(relative_to=__file__)\n\n"
            f"2. Enable parent directory search in your configuration:\n"
            f"   [tool.setuptools_scm]\n"
            f"   search_parent_directories = true\n\n"
            f"3. Change your working directory to the repository root: {scm_parent}\n\n"
            f"4. Set the root explicitly in your configuration:\n"
            f"   [tool.setuptools_scm]\n"
            f'   root = "{scm_parent}"\n\n'
            "For more information, see: https://setuptools-scm.readthedocs.io/en/latest/config/"
        )
    else:
        # No SCM repository found in parent directories either
        error_msg = (
            base_error
            + "Make sure you're either building from a fully intact git repository "
            "or PyPI tarballs. Most other sources (such as GitHub's tarballs, a "
            "git checkout without the .git folder) don't contain the necessary "
            "metadata and will not work.\n\n"
            "For example, if you're using pip, instead of "
            "https://github.com/user/proj/archive/master.zip "
            "use git+https://github.com/user/proj.git#egg=proj\n\n"
            "Alternatively, set the version with the environment variable "
            "SETUPTOOLS_SCM_PRETEND_VERSION_FOR_${NORMALIZED_DIST_NAME} as described "
            "in https://setuptools-scm.readthedocs.io/en/latest/config/"
        )

    raise LookupError(error_msg)


def get_version(
    root: _t.PathT = ".",
    version_scheme: _t.VERSION_SCHEME = _config.DEFAULT_VERSION_SCHEME,
    local_scheme: _t.VERSION_SCHEME = _config.DEFAULT_LOCAL_SCHEME,
    write_to: _t.PathT | None = None,
    write_to_template: str | None = None,
    version_file: _t.PathT | None = None,
    version_file_template: str | None = None,
    relative_to: _t.PathT | None = None,
    tag_regex: str | Pattern[str] = _config.DEFAULT_TAG_REGEX,
    parentdir_prefix_version: str | None = None,
    fallback_version: str | None = None,
    fallback_root: _t.PathT = ".",
    parse: Any | None = None,
    git_describe_command: _t.CMD_TYPE | None = None,
    dist_name: str | None = None,
    version_cls: Any | None = None,
    normalize: bool = True,
    search_parent_directories: bool = False,
    scm: dict[str, Any] | None = None,
) -> str:
    """
    If supplied, relative_to should be a file from which root may
    be resolved. Typically called by a script or module that is not
    in the root of the repository to direct setuptools-scm to the
    root of the repository by supplying ``__file__``.
    """

    version_cls = _validate_version_cls(version_cls, normalize)
    del normalize
    tag_regex = parse_tag_regex(tag_regex)

    # Handle scm parameter by converting it to ScmConfiguration
    if scm is not None:
        scm_config = _config.ScmConfiguration.from_data(scm)
    else:
        scm_config = _config.ScmConfiguration()

    # Remove scm from locals() since we handle it separately
    config_params = locals().copy()
    config_params.pop("scm", None)
    config_params.pop("scm_config", None)

    config = _config.Configuration(scm=scm_config, **config_params)
    maybe_version = _get_version(config, force_write_version_files=True)

    if maybe_version is None:
        _version_missing(config)
    return maybe_version


def parse_tag_regex(tag_regex: str | Pattern[str]) -> Pattern[str]:
    if isinstance(tag_regex, str):
        if tag_regex == "":
            warnings.warn(EMPTY_TAG_REGEX_DEPRECATION)
            return _config.DEFAULT_TAG_REGEX
        else:
            return re.compile(tag_regex)
    else:
        return tag_regex
