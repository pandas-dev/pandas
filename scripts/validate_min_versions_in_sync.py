#!/usr/bin/env python3
"""
Check pandas required and optional dependencies are synced across:

ci/deps/actions-.*-minimum_versions.yaml
pandas/compat/_optional.py
setup.cfg

TODO: doc/source/getting_started/install.rst

This is meant to be run as a pre-commit hook - to run it manually, you can do:

    pre-commit run validate-min-versions-in-sync --all-files
"""
from __future__ import annotations

import os
import pathlib
import sys

import yaml

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

DOC_PATH = pathlib.Path("doc/source/getting_started/install.rst").resolve()
CI_PATH = next(
    pathlib.Path("ci/deps").absolute().glob("actions-*-minimum_versions.yaml")
)
CODE_PATH = pathlib.Path("pandas/compat/_optional.py").resolve()
SETUP_PATH = pathlib.Path("pyproject.toml").resolve()
YAML_PATH = pathlib.Path("ci/deps")
ENV_PATH = pathlib.Path("environment.yml")
EXCLUDE_DEPS = {"tzdata", "blosc"}
# pandas package is not available
# in pre-commit environment
sys.path.append("pandas/compat")
sys.path.append("pandas/util")
import _exceptions
import version

sys.modules["pandas.util.version"] = version
sys.modules["pandas.util._exceptions"] = _exceptions
import _optional


def pin_min_versions_to_ci_deps():
    all_yaml_files = list(YAML_PATH.iterdir())
    all_yaml_files.append(ENV_PATH)
    with open(SETUP_PATH, "rb") as toml_f:
        toml_dic = tomllib.load(toml_f)
    for curr_file in all_yaml_files:
        with open(curr_file) as yaml_f:
            data = update_yaml_file_version(yaml_f, toml_dic)
            os.remove(curr_file)
            with open(curr_file, "w") as f:
                f.write(data)


def update_yaml_file_version(yaml_file, toml_dic):
    exclusion_list = {
        "python=3.8[build=*_pypy]": None,
    }
    toml_deps = {}
    toml_dependencies = set(toml_dic["project"]["optional-dependencies"]["all"])
    for dependency in toml_dependencies:
        toml_package, toml_version = dependency.strip().split(">=")
        toml_deps[toml_package] = toml_version
    data = yaml_file.read()
    yaml_f = yaml.safe_load(data)
    yaml_deps = yaml_f["dependencies"]
    res = []
    [res.append(x) for x in yaml_deps if x not in res]
    yaml_deps = res
    for dep_line in yaml_deps:
        replace_text = operator = yaml_left = ""
        search_text = str(dep_line)
        if str(dep_line) in exclusion_list:
            continue
        if ">=" in dep_line:
            operator = ">="
        elif "==" in dep_line:
            operator = "=="
        elif "=" in dep_line:
            operator = "="
        elif "<" in dep_line:
            operator = "<"
        elif ">" in dep_line:
            operator = ">"
        else:
            operator = ""
        if operator == "":
            yaml_package, yaml_version = str(dep_line).strip(), None
            yaml_left = yaml_package
        else:
            yaml_package, yaml_version = str(dep_line).strip().split(operator)
            if "<" in operator or ">" in operator:
                if yaml_package in toml_deps:
                    if version.parse(yaml_version) <= version.parse(
                        toml_deps[yaml_package]
                    ):
                        if "tzdata" in yaml_package:
                            yaml_left = yaml_package + operator
                        else:
                            yaml_left = yaml_package
                    else:
                        yaml_left = str(dep_line) + ", "
            else:
                yaml_left = yaml_package + operator
        if yaml_package in toml_deps:
            if ">" in yaml_left or "<" in yaml_left:
                if "," in yaml_left:
                    # ex: "numpy<1.24.0," + ">=" + "1.2"
                    replace_text = yaml_left + ">=" + toml_deps[yaml_package]
                else:
                    replace_text = yaml_left + toml_deps[yaml_package]
            # update yaml package version to TOML min version
            elif yaml_version is not None:
                if version.parse(toml_deps[yaml_package]) > version.parse(yaml_version):
                    # ex: "hypothesis>=" + "6.34.2"
                    replace_text = yaml_left + toml_deps[yaml_package]
                elif version.parse(toml_deps[yaml_package]) == version.parse(
                    yaml_version
                ):
                    replace_text = dep_line
            else:
                # ex: "hypothesis + ">=" + 6.34.2"
                replace_text = yaml_package + ">=" + toml_deps[yaml_package]
            data = data.replace(search_text, replace_text)
    return data


def get_versions_from_code() -> dict[str, str]:
    """Min versions for checking within pandas code."""
    install_map = _optional.INSTALL_MAPPING
    versions = _optional.VERSIONS
    for item in EXCLUDE_DEPS:
        versions.pop(item, None)
    return {install_map.get(k, k).casefold(): v for k, v in versions.items()}


def get_versions_from_ci(content: list[str]) -> tuple[dict[str, str], dict[str, str]]:
    """Min versions in CI job for testing all optional dependencies."""
    # Don't parse with pyyaml because it ignores comments we're looking for
    seen_required = False
    seen_optional = False
    seen_test = False
    required_deps = {}
    optional_deps = {}
    for line in content:
        if "# test dependencies" in line:
            seen_test = True
        elif seen_test and "- pytest>=" in line:
            # Only grab pytest
            package, version = line.strip().split(">=")
            package = package[2:]
            optional_deps[package.casefold()] = version
        elif "# required dependencies" in line:
            seen_required = True
        elif "# optional dependencies" in line:
            seen_optional = True
        elif "- pip:" in line:
            continue
        elif seen_required and line.strip():
            if "==" in line:
                package, version = line.strip().split("==")

            else:
                package, version = line.strip().split("=")
            package = package[2:]
            if package in EXCLUDE_DEPS:
                continue
            if not seen_optional:
                required_deps[package.casefold()] = version
            else:
                optional_deps[package.casefold()] = version
    return required_deps, optional_deps


def get_versions_from_toml() -> dict[str, str]:
    """Min versions in pyproject.toml for pip install pandas[extra]."""
    install_map = _optional.INSTALL_MAPPING
    optional_dependencies = {}
    with open(SETUP_PATH, "rb") as pyproject_f:
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


def main():
    pin_min_versions_to_ci_deps()
    with open(CI_PATH, encoding="utf-8") as f:
        _, ci_optional = get_versions_from_ci(f.readlines())
    code_optional = get_versions_from_code()
    setup_optional = get_versions_from_toml()

    diff = (ci_optional.items() | code_optional.items() | setup_optional.items()) - (
        ci_optional.items() & code_optional.items() & setup_optional.items()
    )

    if diff:
        packages = {package for package, _ in diff}
        out = sys.stdout
        out.write(
            f"The follow minimum version differences were found between  "
            f"{CI_PATH}, {CODE_PATH} AND {SETUP_PATH}. "
            f"Please ensure these are aligned: \n\n"
        )

        for package in packages:
            out.write(
                f"{package}\n"
                f"{CI_PATH}: {ci_optional.get(package, 'Not specified')}\n"
                f"{CODE_PATH}: {code_optional.get(package, 'Not specified')}\n"
                f"{SETUP_PATH}: {setup_optional.get(package, 'Not specified')}\n\n"
            )
        sys.exit(1)
    sys.exit(0)


if __name__ == "__main__":
    main()
