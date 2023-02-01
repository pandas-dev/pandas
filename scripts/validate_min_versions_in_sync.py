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
    toml_dependencies = get_versions_from_toml()
    print("TOML Dependencies", toml_dependencies)
    # if toml package is in yaml file, then run condition check
    # TODO: Add environment.yml to "all_yaml_files"
    all_yaml_files = list(YAML_PATH.iterdir())
    yaml_dependencies = {}
    # iterate through each yaml file
    for curr_file in all_yaml_files:
        with open(curr_file, "rb") as yaml_f:
            # all stuff including dependencies
            yaml_file = yaml.safe_load(yaml_f)
            print("\nCurrent YAML File Path:", curr_file)
            # just yaml deps
            yaml_deps = yaml_file["dependencies"]
            print("YAML Dictionary:", yaml_deps)
            # iterate through package/version dependency string
            for dependency in yaml_deps:
                # clean the string. extract dictionary data: yaml package, version
                if ">=" in dependency:
                    yaml_package, yaml_version = str(dependency).strip().split(">=")
                elif "=" in dependency:
                    package, version = str(dependency).strip().split("=")
                else:
                    package, version = str(dependency), "None"
                # comparison between YAML/TOML
                if ">=" in dependency or "=" in dependency:
                    if yaml_package in toml_dependencies:
                        if toml_dependencies[yaml_package] > yaml_version:
                            # update yaml package version to toml min version
                            pass
                else:
                    if yaml_package in toml_dependencies:
                        # update yaml package version to toml min version
                        # use ">=" operator
                        pass

                # put extracted data into yaml dictionary
                yaml_dependencies[package] = version
            print()
            myKeys = list(yaml_dependencies.keys())
            myKeys.sort()
            sorted_dict = {i: yaml_dependencies[i] for i in myKeys}

            print("Sorted YAML Dependencies:")
            for item in sorted_dict.items():
                print(item)
            print("===================")
    return yaml_dependencies


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
    print()
    myKeys = list(optional_dependencies.keys())
    myKeys.sort()
    sorted_dict = {i: optional_dependencies[i] for i in myKeys}
    print("Sorted TOML Dependencies:")
    for item in sorted_dict.items():
        print(item)

    return optional_dependencies


def main():
    yaml_dependencies = pin_min_versions_to_ci_deps()
    print("YAML Dependencies:", yaml_dependencies)

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
