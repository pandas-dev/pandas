#!/usr/bin/env python3
"""
Check that every dependency the test suite gates on is installed somewhere in CI.

``pytest.importorskip("foo")`` and ``td.skip_if_no("foo")`` silently skip when
``foo`` is missing.  If ``foo`` is not installed in *any* environment that runs
the test suite, the tests behind that gate never run anywhere and the skip is
invisible -- locally they pass, because a development environment installs more
than CI does.

This is meant to be run as a pre-commit hook - to run it manually, you can do:

    pre-commit run validate-test-dependencies --all-files
"""

from __future__ import annotations

import argparse
import pathlib
import re
import sys
import tomllib
from typing import TYPE_CHECKING

import yaml

if TYPE_CHECKING:
    from collections.abc import (
        Iterable,
        Sequence,
    )

BASE_PATH = pathlib.Path(__file__).parents[1]
PIXI_PATH = BASE_PATH / "pixi.toml"
WORKFLOW_PATH = BASE_PATH / ".github/workflows/unit-tests.yml"
PANDAS_PATH = BASE_PATH / "pandas"

# skip_if_no is defined here, and its docstring shows the usage. Scanning it
# would pick up the placeholder name from that example.
DEFINITION_SITE = BASE_PATH / "pandas/util/_test_decorators.py"

GATE_PATTERN = re.compile(
    r"""importorskip\(\s*["']([A-Za-z0-9_.]+)["']"""
    r"""|skip_if_no\(\s*["']([A-Za-z0-9_.]+)["']"""
)

# Import name -> conda package name, for the cases where they differ.
IMPORT_TO_CONDA = {
    "IPython": "ipython",
    "adbc_driver_postgresql": "adbc-driver-postgresql",
    "adbc_driver_sqlite": "adbc-driver-sqlite",
    "bs4": "beautifulsoup4",
    "odf": "odfpy",
    "pandas_datareader": "pandas-datareader",
    "python_calamine": "python-calamine",
    "sklearn": "scikit-learn",
    "tables": "pytables",
    "yaml": "pyyaml",
}

# Modules that are legitimately absent from pixi.toml, with the reason why.
ALLOWED_UNDECLARED = {
    "moto": "CI serves moto from a container, see PANDAS_MOTO_URL",
}

MESSAGE = (
    "{module!r} gates tests with importorskip/skip_if_no, but {package!r} is "
    "not installed in any environment that runs the test suite, so these "
    "tests never run:"
)


def environments_running_tests() -> set[str]:
    """
    Names of the pixi environments that the unit test workflow runs pytest in.
    """
    with open(WORKFLOW_PATH, encoding="utf-8") as file_handle:
        workflow = yaml.safe_load(file_handle)

    environments = set()
    for job in workflow["jobs"].values():
        matrix = job.get("strategy", {}).get("matrix", {})
        environments.update(matrix.get("environment", []))
        for entry in matrix.get("include", []):
            if "environment" in entry:
                environments.add(entry["environment"])
    return environments


def declared_packages(pixi: dict, environments: Iterable[str]) -> set[str]:
    """
    Every package installed in at least one of ``environments``.
    """
    packages = set(pixi.get("dependencies", {}))
    for environment in environments:
        specification = pixi["environments"][environment]
        if isinstance(specification, dict):
            features = specification["features"]
        else:
            features = specification
        for feature in features:
            block = pixi["feature"].get(feature, {})
            packages.update(block.get("dependencies", {}))
            packages.update(block.get("pypi-dependencies", {}))
    return packages


def gated_modules() -> dict[str, set[str]]:
    """
    Map each gated top level module to the files that gate on it.
    """
    modules: dict[str, set[str]] = {}
    for path in sorted(PANDAS_PATH.rglob("*.py")):
        if path == DEFINITION_SITE:
            continue
        for match in GATE_PATTERN.finditer(path.read_text(encoding="utf-8")):
            module = (match.group(1) or match.group(2)).split(".")[0]
            modules.setdefault(module, set()).add(str(path.relative_to(BASE_PATH)))
    return modules


def validate_test_dependencies() -> int:
    with open(PIXI_PATH, "rb") as file_handle:
        pixi = tomllib.load(file_handle)
    declared = declared_packages(pixi, environments_running_tests())

    ret = 0
    for module, paths in sorted(gated_modules().items()):
        if module in ALLOWED_UNDECLARED or module in sys.stdlib_module_names:
            continue
        package = IMPORT_TO_CONDA.get(module, module)
        if package in declared:
            continue
        print(MESSAGE.format(module=module, package=package))
        for path in sorted(paths):
            print(f"    {path}")
        ret = 1
    return ret


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args(argv)
    return validate_test_dependencies()


if __name__ == "__main__":
    sys.exit(main())
