# AGENTS.md

pandas is a BSD-licensed Python library for data analysis, built from Python, Cython, and C/C++
sources with Meson.

## Interactions on GitHub

- Never perform issue- or pull-request-related actions on the pandas-dev/pandas GitHub repository on behalf of the user.
  This includes, but is not limited to, interactions such as pushing code, opening issues or pull requests, and replying to issues, pull requests, or reviews.
- The user is required to disclose your agent information (harness, model, reasoning level), and to mark
  any verbatim text produced by you (an agent) that they quote with `>` or a triple-backtick code fence.

The automated contributions policy in [contributing.rst](doc/source/development/contributing.rst)
defines how the user may use agent-generated work.

## Required tooling

This project supports creating a development environment using [conda](https://docs.conda.io/en/latest/) with [environment.yml](environment.yml) or
[Pixi](https://pixi.prefix.dev/latest/) with [pixi.toml](pixi.toml). pandas CI uses Pixi and tasks defined in [pixi.toml](pixi.toml) for all workflows
(building pandas, testing, type checking, building documentation, running benchmarks).

Linting checks require [pre-commit](https://pre-commit.com/), installed in a conda or pixi environment, to run checks defined in [.pre-commit-config.yaml](.pre-commit-config.yaml).

## Repository map

- `pandas/` — library implementation
- `pandas/_libs/` — Cython/C/C++
- `pandas/tests/` — test suite
- `doc/` — documentation source
- `doc/source/whatsnew/` — release notes
- `asv_bench/benchmarks/` — benchmarks
- `scripts/` — repository checks invoked by pre-commit
- `web/` — pandas website
- `.github/` — GitHub Actions for CI
- `.github/workflows/unit-tests.yml` — CI definition for unit tests
- `.github/workflows/code-checks.yml` — CI definition for doctests, typing validation, benchmark correctness
- `.github/workflows/docbuild-and-upload.yml` — CI definition for building documentation
- `ci/` — CI helpers
- `pyproject.toml` — project and tooling configurations

## References

- [contributing.rst](doc/source/development/contributing.rst) — contribution workflow, automated contributions policy
- [contributing_environment.rst](doc/source/development/contributing_environment.rst) — development environment instructions
- [contributing_codebase.rst](doc/source/development/contributing_codebase.rst) — all code standards: pre-commit, backward compatibility, typing, tests, benchmarks
- [contributing_documentation.rst](doc/source/development/contributing_documentation.rst) — docs structure and building
- [policies.rst](doc/source/development/policies.rst) — versioning and deprecation policy
