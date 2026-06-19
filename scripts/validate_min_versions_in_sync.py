#!/usr/bin/env python3
"""
Check pandas required and optional dependencies are synced across:

pixi.toml
pandas/compat/_optional.py
pyproject.toml
environment.yml

TODO: doc/source/getting_started/install.rst

This is meant to be run as a pre-commit hook - to run it manually, you can do:

    pre-commit run validate-min-versions-in-sync --all-files
"""

from __future__ import annotations

import pathlib
import sys
import tomllib
from typing import Any

import yaml

from scripts.generate_pip_deps_from_conda import CONDA_TO_PIP

BASE_PATH = pathlib.Path(__file__).parents[1]
DOC_PATH = (BASE_PATH / "doc/source/getting_started/install.rst").resolve()
PIXI_PATH = (BASE_PATH / "pixi.toml").resolve()
CODE_PATH = (BASE_PATH / "pandas/compat/_optional.py").resolve()
PYPROJECT_PATH = (BASE_PATH / "pyproject.toml").resolve()
ENV_PATH = BASE_PATH / "environment.yml"
EXCLUDE_DEPS = {"blosc", "tzdata", "pyqt", "pyqt5"}
ADDITIONAL_PIXI_OPTIONAL_DEPENDENCIES = {
    ("feature", "pyarrow", "dependencies"): None,
    ("feature", "test-base", "dependencies"): {"hypothesis", "pytest"},
    ("feature", "test-clipboard", "dependencies"): {"qtpy"},
}
# pandas package is not available
# in pre-commit environment
sys.path.append(str(BASE_PATH / "pandas/compat"))
sys.path.append(str(BASE_PATH / "pandas/util"))
import _exceptions
import version

sys.modules["pandas.util.version"] = version
sys.modules["pandas.util._exceptions"] = _exceptions
import _optional


def pin_min_versions_to_environment_yml() -> int:
    """
    Pin minimum versions to environment.yml dependencies.

    Pip dependencies are not pinned.
    """
    toml_dependencies = {}
    with open(PYPROJECT_PATH, "rb") as toml_f:
        toml_dependencies = tomllib.load(toml_f)
    ret = 0
    with open(ENV_PATH, encoding="utf-8") as yaml_f:
        yaml_start_data = yaml_f.read()
    yaml_file = yaml.safe_load(yaml_start_data)
    yaml_dependencies = yaml_file["dependencies"]
    yaml_map = get_yaml_map_from(yaml_dependencies)
    toml_map = get_toml_map_from(toml_dependencies)
    yaml_result_data = pin_min_versions_to_yaml_file(
        yaml_map, toml_map, yaml_start_data
    )
    if yaml_result_data != yaml_start_data:
        with open(ENV_PATH, "w", encoding="utf-8") as f:
            f.write(yaml_result_data)
            ret |= 1
    return ret


def get_toml_map_from(toml_dic: dict[str, Any]) -> dict[str, str]:
    toml_deps = {}
    toml_dependencies = set(toml_dic["project"]["optional-dependencies"]["all"])
    for dependency in toml_dependencies:
        toml_package, toml_version = dependency.strip().split(">=")
        toml_deps[toml_package] = toml_version
    return toml_deps


def get_operator_from(dependency: str) -> str | None:
    if "<=" in dependency:
        operator = "<="
    elif ">=" in dependency:
        operator = ">="
    elif "=" in dependency:
        operator = "="
    elif ">" in dependency:
        operator = ">"
    elif "<" in dependency:
        operator = "<"
    else:
        operator = None
    return operator


def get_yaml_map_from(
    yaml_dic: list[str | dict[str, list[str]]],
) -> dict[str, list[str] | None]:
    yaml_map: dict[str, list[str] | None] = {}
    for dependency in yaml_dic:
        if isinstance(dependency, dict) or dependency in yaml_map:
            continue
        search_text = str(dependency)
        operator = get_operator_from(search_text)
        if "," in dependency:
            yaml_dependency, yaml_version1 = search_text.split(",")
            operator = get_operator_from(yaml_dependency)
            assert operator is not None
            yaml_package, yaml_version2 = yaml_dependency.split(operator)
            yaml_version2 = operator + yaml_version2
            yaml_map[yaml_package] = [yaml_version1, yaml_version2]
        elif operator is not None:
            yaml_package, yaml_version = search_text.split(operator)
            yaml_version = operator + yaml_version
            yaml_map[yaml_package] = [yaml_version]
        else:
            yaml_package, yaml_version = search_text.strip(), None
            yaml_map[yaml_package] = yaml_version
    return yaml_map


def clean_version_list(
    yaml_versions: list[str], toml_version: version.Version
) -> list[str]:
    for i in range(len(yaml_versions)):
        yaml_version = yaml_versions[i].strip()
        operator = get_operator_from(yaml_version)
        assert operator is not None
        if "<=" in operator or ">=" in operator:
            yaml_version = yaml_version[2:]
        else:
            yaml_version = yaml_version[1:]
        yaml_version = version.parse(yaml_version)
        if yaml_version < toml_version:
            yaml_versions[i] = "-" + str(yaml_version)
        elif yaml_version >= toml_version:
            if ">" in operator:
                yaml_versions[i] = "-" + str(yaml_version)
    return yaml_versions


def pin_min_versions_to_yaml_file(
    yaml_map: dict[str, list[str] | None], toml_map: dict[str, str], yaml_file_data: str
) -> str:
    data = yaml_file_data
    for yaml_package, yaml_versions in yaml_map.items():
        old_dep = yaml_package
        if yaml_versions is not None:
            old_dep = old_dep + ", ".join(yaml_versions)
        if CONDA_TO_PIP.get(yaml_package, yaml_package) in toml_map:
            min_dep = toml_map[CONDA_TO_PIP.get(yaml_package, yaml_package)]
        elif yaml_package in toml_map:
            min_dep = toml_map[yaml_package]
        else:
            continue
        if yaml_versions is None:
            new_dep = old_dep + ">=" + min_dep
            data = data.replace(old_dep, new_dep, 1)
            continue
        toml_version = version.parse(min_dep)
        yaml_versions_list = clean_version_list(yaml_versions, toml_version)
        cleaned_yaml_versions = [x for x in yaml_versions_list if "-" not in x]
        new_dep = yaml_package
        for clean_yaml_version in cleaned_yaml_versions:
            new_dep += clean_yaml_version + ", "
        operator = get_operator_from(new_dep)
        if operator != "=":
            new_dep += ">=" + min_dep
        else:
            new_dep = new_dep[:-2]
        data = data.replace(old_dep, new_dep)
    return data


def get_versions_from_code() -> dict[str, str]:
    """Min versions for checking within pandas code."""
    install_map = _optional.INSTALL_MAPPING
    inverse_install_map = {v: k for k, v in install_map.items()}
    versions = _optional.VERSIONS
    for item in EXCLUDE_DEPS:
        item = inverse_install_map.get(item, item)
        versions.pop(item, None)
    return {install_map.get(k, k).casefold(): v for k, v in versions.items()}


def get_versions_from_ci(content: list[str]) -> tuple[dict[str, str], dict[str, str]]:
    """Min versions in pixi.toml for testing all optional dependencies."""
    separator = "" if any(line.endswith(("\n", "\r")) for line in content) else "\n"
    pixi_toml = tomllib.loads(separator.join(content))
    install_map = _optional.INSTALL_MAPPING
    required_deps = _get_required_dependencies_from_pixi_content(content, pixi_toml)
    optional_deps = _get_dependencies_from_pixi_table(
        pixi_toml["feature"]["optional-dependencies"]["dependencies"],
        install_map,
    )
    required_deps.update(
        _get_dependencies_from_pixi_table(
            pixi_toml["feature"]["numpy"]["dependencies"],
            install_map,
        )
    )
    for keys, packages in ADDITIONAL_PIXI_OPTIONAL_DEPENDENCIES.items():
        dependencies = pixi_toml
        for key in keys:
            dependencies = dependencies[key]
        optional_deps.update(
            _get_dependencies_from_pixi_table(dependencies, install_map, packages)
        )
    return required_deps, optional_deps


def _get_required_dependencies_from_pixi_content(
    content: list[str], pixi_toml: dict[str, Any]
) -> dict[str, str]:
    # pixi.toml keeps build and required deps in the same TOML table.
    required_deps = {}
    in_dependencies = False
    in_required_dependencies = False
    for line in content:
        line = line.strip()
        if not line:
            continue
        if line.startswith("["):
            in_dependencies = line == "[dependencies]"
            in_required_dependencies = False
            continue
        if not in_dependencies:
            continue
        if line.casefold() == "# required dependencies":
            in_required_dependencies = True
            continue
        if line.startswith("#"):
            in_required_dependencies = False
            continue
        if not in_required_dependencies:
            continue
        package = line.split("=", maxsplit=1)[0].strip()
        version = _get_version_from_pixi_spec(pixi_toml["dependencies"][package])
        if version is not None and package not in EXCLUDE_DEPS:
            required_deps[package.casefold()] = version
    return required_deps


def _get_dependencies_from_pixi_table(
    dependencies: dict[str, str | dict[str, Any]],
    install_map: dict[str, str],
    packages: set[str] | None = None,
) -> dict[str, str]:
    versions = {}
    for package, spec in dependencies.items():
        package = install_map.get(package, package).casefold()
        if package in EXCLUDE_DEPS or (
            packages is not None and package not in packages
        ):
            continue
        version = _get_version_from_pixi_spec(spec)
        if version is not None:
            versions[package] = version
    return versions


def _get_version_from_pixi_spec(spec: str | dict[str, Any]) -> str | None:
    if isinstance(spec, dict):
        spec = spec.get("version", "")
    if spec in {"", "*"}:
        return None
    for operator in ("==", ">=", "="):
        if spec.startswith(operator):
            version = spec.removeprefix(operator).split(",", maxsplit=1)[0]
            if "*" not in version:
                return version
    return None


def get_versions_from_toml() -> dict[str, str]:
    """Min versions in pyproject.toml for pip install pandas[extra]."""
    install_map = _optional.INSTALL_MAPPING
    optional_dependencies = {}
    with open(PYPROJECT_PATH, "rb") as pyproject_f:
        pyproject_toml = tomllib.load(pyproject_f)
        opt_deps = pyproject_toml["project"]["optional-dependencies"]
        dependencies = set(opt_deps["all"])

        # remove pytest plugin dependencies
        pytest_plugins = {dep for dep in opt_deps["test"] if dep.startswith("pytest-")}
        dependencies = dependencies.difference(pytest_plugins)

    for dependency in dependencies:
        package, version = dependency.strip().split(">=")
        optional_dependencies[install_map.get(package, package).casefold()] = version

    for item in EXCLUDE_DEPS:
        optional_dependencies.pop(item, None)
    return optional_dependencies


def main() -> int:
    ret = 0
    ret |= pin_min_versions_to_environment_yml()
    with open(PIXI_PATH, encoding="utf-8") as f:
        _, ci_optional = get_versions_from_ci(f.readlines())
    code_optional = get_versions_from_code()
    pyproject_optional = get_versions_from_toml()

    diff = (
        ci_optional.items() | code_optional.items() | pyproject_optional.items()
    ) - (ci_optional.items() & code_optional.items() & pyproject_optional.items())

    if diff:
        packages = {package for package, _ in diff}
        out = sys.stdout
        out.write(
            f"The follow minimum version differences were found between  "
            f"{PIXI_PATH}, {CODE_PATH} AND {PYPROJECT_PATH}. "
            f"Please ensure these are aligned: \n\n"
        )

        for package in packages:
            out.write(
                f"{package}\n"
                f"{PIXI_PATH}: {ci_optional.get(package, 'Not specified')}\n"
                f"{CODE_PATH}: {code_optional.get(package, 'Not specified')}\n"
                f"{PYPROJECT_PATH}: {pyproject_optional.get(package, 'Not specified')}\n\n"  # noqa: E501
            )
        ret |= 1
    return ret


if __name__ == "__main__":
    sys.exit(main())
