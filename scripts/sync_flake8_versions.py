"""
Check that the flake8 pins are the same in:

- environment.yml
- .pre-commit-config.yaml, in the flake8 hook
- .pre-commit-config.yaml, in the additional dependencies of the yesqa hook

The flake8 hook revision in .pre-commit-config.yaml is taken as the reference revision.

Usage: either

- ``python scripts/sync_flake8_versions.py``, or
- ``pre-commit run sync-flake8-versions --all-files``.
"""
import sys
from typing import (
    Any,
    Mapping,
    NamedTuple,
    Sequence,
    Tuple,
    TypeVar,
)

import yaml


class Revisions(NamedTuple):
    precommit_rev: str
    precommit_yesqa_rev: str
    environment_rev: str


YamlMapping = Mapping[str, Any]
Repo = TypeVar("Repo", bound=YamlMapping)


def _get_repo_hook(repos: Sequence[Repo], hook_name: str) -> Tuple[Repo, YamlMapping]:
    for repo in repos:
        for hook in repo["hooks"]:
            if hook["id"] == hook_name:
                return repo, hook
    else:
        raise RuntimeError(f"Repo with hook {hook_name} not found")


def get_revisions(precommit_config: YamlMapping, environment: YamlMapping) -> Revisions:
    repos = precommit_config["repos"]
    flake8_repo, _ = _get_repo_hook(repos, "flake8")
    precommit_rev = flake8_repo["rev"]

    _, yesqa_hook = _get_repo_hook(repos, "yesqa")
    additional_dependencies = yesqa_hook.get("additional_dependencies", [])
    for dep in additional_dependencies:
        if "==" in dep:
            pkg, rev = dep.split("==", maxsplit=1)
            if pkg == "flake8":
                precommit_yesqa_rev = rev
                break
    else:
        raise RuntimeError(
            "flake8 not found, or not pinned, in additional dependencies of yesqa "
            "hook in .pre-commit-config.yaml"
        )

    deps = environment["dependencies"]
    for dep in deps:
        if isinstance(dep, str) and "=" in dep:
            pkg, rev = dep.split("=", maxsplit=1)
            if pkg == "flake8":
                environment_rev = rev
                break
    else:
        raise RuntimeError("flake8 not found, or not pinned, in environment.yml")

    return Revisions(precommit_rev, precommit_yesqa_rev, environment_rev)


if __name__ == "__main__":
    with open(".pre-commit-config.yaml") as fd:
        precommit_config = yaml.safe_load(fd)
    with open("environment.yml") as fd:
        environment = yaml.safe_load(fd)

    revisions = get_revisions(precommit_config, environment)

    if revisions.environment_rev != revisions.precommit_rev:
        sys.stdout.write(
            f"flake8 pin in environment.yml is {revisions.environment_rev}, "
            f"should be {revisions.precommit_rev}\n"
        )
        sys.exit(1)

    if revisions.precommit_yesqa_rev != revisions.precommit_rev:
        sys.stdout.write(
            f"flake8 pin in yesqa is {revisions.precommit_yesqa_rev}, "
            f"should be {revisions.precommit_rev}\n"
        )
        sys.exit(1)

    sys.exit(0)
