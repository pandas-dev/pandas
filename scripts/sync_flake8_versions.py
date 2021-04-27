"""
Check that the flake8 (and pandas-dev-flaker) pins are the same in:

- environment.yml
- .pre-commit-config.yaml, in the flake8 hook
- .pre-commit-config.yaml, in the additional dependencies of the yesqa hook

The flake8 hook revision in .pre-commit-config.yaml is taken as the reference revision.

Usage: either

- ``python scripts/sync_flake8_versions.py``, or
- ``pre-commit run sync-flake8-versions --all-files``.
"""
from __future__ import annotations

from dataclasses import (
    dataclass,
    replace,
)
import sys
from typing import (
    Any,
    Mapping,
    Sequence,
    TypeVar,
)

import yaml


@dataclass
class Revision:
    name: str
    compare: str
    version: str


@dataclass
class Revisions:
    name: str
    pre_commit: Revision | None = None
    yesqa: Revision | None = None
    environment: Revision | None = None


YamlMapping = Mapping[str, Any]
Repo = TypeVar("Repo", bound=YamlMapping)

COMPARE = ("<=", "==", ">=", "<", ">", "=")


def _get_repo_hook(repos: Sequence[Repo], hook_name: str) -> tuple[Repo, YamlMapping]:
    for repo in repos:
        for hook in repo["hooks"]:
            if hook["id"] == hook_name:
                return repo, hook
    else:  # pragma: no cover
        raise RuntimeError(f"Repo with hook {hook_name} not found")


def _conda_to_pip_compat(dep):
    if dep.compare == "=":
        return replace(dep, compare="==")
    else:
        return dep


def _validate_additional_dependencies(
    flake8_additional_dependencies,
    yesqa_additional_dependencies,
    environment_additional_dependencies,
) -> None:
    for dep in flake8_additional_dependencies:
        if dep not in yesqa_additional_dependencies:
            sys.stdout.write(
                f"Mismatch of '{dep.name}' version between 'flake8' "
                "and 'yesqa' in '.pre-commit-config.yaml'\n"
            )
            sys.exit(1)
        if dep not in environment_additional_dependencies:
            sys.stdout.write(
                f"Mismatch of '{dep.name}' version between 'enviroment.yml' "
                "and additional dependencies of 'flake8' in '.pre-commit-config.yaml'\n"
            )
            sys.exit(1)


def _validate_revisions(revisions):
    if revisions.environment != revisions.pre_commit:
        sys.stdout.write(
            f"{revisions.name} in 'environment.yml' does not "
            "match in 'flake8' from 'pre-commit'\n"
        )
        sys.exit(1)

    if revisions.yesqa != revisions.pre_commit:
        sys.stdout.write(
            f"{revisions.name} in 'yesqa' does not match "
            "in 'flake8' from 'pre-commit'\n"
        )
        sys.exit(1)


def _process_dependencies(deps):
    for dep in deps:
        if isinstance(dep, str):
            for compare in COMPARE:
                if compare in dep:
                    pkg, rev = dep.split(compare, maxsplit=1)
                    yield _conda_to_pip_compat(Revision(pkg, compare, rev))
                    break
        else:
            yield from _process_dependencies(dep["pip"])


def get_revisions(
    precommit_config: YamlMapping, environment: YamlMapping
) -> tuple[Revisions, Revisions]:
    flake8_revisions = Revisions(name="flake8")
    pandas_dev_flaker_revisions = Revisions(name="pandas-dev-flaker")

    repos = precommit_config["repos"]
    flake8_repo, flake8_hook = _get_repo_hook(repos, "flake8")
    flake8_revisions.pre_commit = Revision("flake8", "==", flake8_repo["rev"])
    flake8_additional_dependencies = []
    for dep in _process_dependencies(flake8_hook.get("additional_dependencies", [])):
        if dep.name == "pandas-dev-flaker":
            pandas_dev_flaker_revisions.pre_commit = dep
        else:
            flake8_additional_dependencies.append(dep)

    _, yesqa_hook = _get_repo_hook(repos, "yesqa")
    yesqa_additional_dependencies = []
    for dep in _process_dependencies(yesqa_hook.get("additional_dependencies", [])):
        if dep.name == "flake8":
            flake8_revisions.yesqa = dep
        elif dep.name == "pandas-dev-flaker":
            pandas_dev_flaker_revisions.yesqa = dep
        else:
            yesqa_additional_dependencies.append(dep)

    environment_dependencies = environment["dependencies"]
    environment_additional_dependencies = []
    for dep in _process_dependencies(environment_dependencies):
        if dep.name == "flake8":
            flake8_revisions.environment = dep
        elif dep.name == "pandas-dev-flaker":
            pandas_dev_flaker_revisions.environment = dep
        else:
            environment_additional_dependencies.append(dep)

    _validate_additional_dependencies(
        flake8_additional_dependencies,
        yesqa_additional_dependencies,
        environment_additional_dependencies,
    )

    for revisions in flake8_revisions, pandas_dev_flaker_revisions:
        _validate_revisions(revisions)


if __name__ == "__main__":
    with open(".pre-commit-config.yaml") as fd:
        precommit_config = yaml.safe_load(fd)
    with open("environment.yml") as fd:
        environment = yaml.safe_load(fd)
    get_revisions(precommit_config, environment)
    sys.exit(0)
