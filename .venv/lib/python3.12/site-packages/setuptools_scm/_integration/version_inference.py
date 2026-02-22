from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING
from typing import Any
from typing import Union

from setuptools import Distribution

from .. import _log

if TYPE_CHECKING:
    from .pyproject_reading import PyProjectData

log = _log.log.getChild("version_inference")


@dataclass
class VersionInferenceConfig:
    """Configuration for version inference."""

    dist_name: str | None
    pyproject_data: PyProjectData | None
    overrides: dict[str, Any] | None

    def apply(self, dist: Distribution) -> None:
        """Apply version inference to the distribution."""
        version_string = infer_version_string(
            self.dist_name,
            self.pyproject_data,  # type: ignore[arg-type]
            self.overrides,
            force_write_version_files=True,
        )
        dist.metadata.version = version_string

        # Mark that this version was set by infer_version if overrides is None (infer_version context)
        if self.overrides is None:
            dist._setuptools_scm_version_set_by_infer = True  # type: ignore[attr-defined]


@dataclass
class VersionInferenceWarning:
    """Error message for user."""

    message: str

    def apply(self, dist: Distribution) -> None:
        """Apply error handling to the distribution."""
        import warnings

        warnings.warn(self.message)


@dataclass(frozen=True)
class VersionInferenceNoOp:
    """No operation result - silent skip."""

    def apply(self, dist: Distribution) -> None:
        """Apply no-op to the distribution."""


VersionInferenceResult = Union[
    VersionInferenceConfig,  # Proceed with inference
    VersionInferenceWarning,  # Show warning
    VersionInferenceNoOp,  # Don't infer (silent)
]


def infer_version_string(
    dist_name: str | None,
    pyproject_data: PyProjectData,
    overrides: dict[str, Any] | None = None,
    *,
    force_write_version_files: bool = False,
) -> str:
    """
    Compute the inferred version string from the given inputs without requiring a
    setuptools Distribution instance. This is a pure helper that simplifies
    integration tests by avoiding file I/O and side effects on a Distribution.

    Parameters:
        dist_name: Optional distribution name (used for overrides and env scoping)
        pyproject_data: Parsed PyProjectData (may be constructed via for_testing())
        overrides: Optional override configuration (same keys as [tool.setuptools_scm])
        force_write_version_files: When True, apply write_to/version_file effects

    Returns:
        The computed version string.
    """
    from .. import _config as _config_module
    from .._get_version_impl import _get_version
    from .._get_version_impl import _version_missing

    config = _config_module.Configuration.from_file(
        dist_name=dist_name, pyproject_data=pyproject_data, **(overrides or {})
    )

    maybe_version = _get_version(
        config, force_write_version_files=force_write_version_files
    )
    if maybe_version is None:
        _version_missing(config)
    return maybe_version


def get_version_inference_config(
    dist_name: str | None,
    current_version: str | None,
    pyproject_data: PyProjectData,
    overrides: dict[str, Any] | None = None,
) -> VersionInferenceResult:
    """
    Determine whether and how to perform version inference.

    Args:
        dist_name: The distribution name
        current_version: Current version if any
        pyproject_data: PyProjectData from parser (None if file doesn't exist)
        overrides: Override configuration (None for no overrides)

    Returns:
        VersionInferenceResult with the decision and configuration
    """

    config = VersionInferenceConfig(
        dist_name=dist_name,
        pyproject_data=pyproject_data,
        overrides=overrides,
    )

    inference_implied = pyproject_data.should_infer() or overrides is not None

    if inference_implied:
        if current_version is None:
            return config
        else:
            return VersionInferenceWarning(
                f"version of {dist_name} already set",
            )
    else:
        return VersionInferenceNoOp()
