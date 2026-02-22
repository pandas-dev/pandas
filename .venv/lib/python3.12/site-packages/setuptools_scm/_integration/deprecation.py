import warnings

from pathlib import Path


def warn_dynamic_version(path: Path, section: str, expression: str) -> None:
    warnings.warn(
        f"{path}: at [{section}]\n"
        f"{expression} is forcing setuptools to override the version setuptools-scm did already set\n"
        "When using setuptools-scm it's invalid to use setuptools dynamic version as well, please remove it.\n"
        "Setuptools-scm is responsible for setting the version, forcing setuptools to override creates errors."
    )


def warn_pyproject_setuptools_dynamic_version(path: Path) -> None:
    warn_dynamic_version(path, "tool.setuptools.dynamic", "version = {attr = ...}")


def warn_setup_cfg_dynamic_version(path: Path) -> None:
    warn_dynamic_version(path, "metadata", "version = attr: ...")
