# AGENTS.md

pandas is a BSD-licensed Python library for data analysis, built from Python, Cython, and C/C++
sources with Meson.

## Interactions on GitHub

- Never perform issue- or pull-request-related actions on the pandas-dev/pandas GitHub repository on behalf of the user.
  This includes, but is not limited to, interactions such as pushing code, opening issues or pull requests, and replying to issues, pull requests, or reviews.
- The user is required to disclose your agent information (harness, model, reasoning level), and to mark
  any verbatim text produced by you (an agent) that they quote with `>` or a triple-backtick code fence.

The automated contributions policy in
[contributing.rst](doc/source/development/contributing.rst) defines how the user may use agent-generated work.

## Required tooling

[pixi](https://pixi.prefix.dev/latest/) or [conda](https://docs.conda.io/en/latest/) are used to manage a development environment and its dependencies.

| Stack | Setup | Running commands |
| --- | --- | --- |
| conda (recommended in the contributor docs) | `conda env create --file environment.yml` | `conda run -n pandas-dev <command>` |
| pixi (used by CI, lockfile committed) | `pixi install --environment <env>` | `pixi run --environment <env> <command>` |

For conda, use a pre-existing `pandas-dev` environment if available (discovered by running `conda env list`).
The conda environment contains dependencies to run all applicable tasks.

Relevant pixi environments, defined in [pixi.toml](pixi.toml):

| Environment | Use |
| --- | --- |
| `py311` | default for building and running tests locally |
| `typing` | `pre-commit`, `mypy`, `pyright`, `stubtest` |
| `doctests` | doctests |
| `docs-and-web` | Sphinx docs, whatsnew, website |
| `asv` | benchmarks |

The pixi environments contain dependencies to build pandas and run an applicable task. Pivot to a different pixi environment
to run another, distinct task.

## Task and validation navigation

Pick the row matching the work in progress and read the linked guide if needed. The pixi commands mirror what CI runs;
the conda commands are scoped more narrowly to the files you changed.

| Task | Run when | Documentation | conda command | pixi command |
| --- | --- | --- | --- | --- |
| Build / rebuild | Before any other task, and again after editing Cython, C, or C++ sources | [contributing_environment.rst](doc/source/development/contributing_environment.rst) | `conda run -n pandas-dev python -m pip install --verbose --editable . --no-build-isolation` | `pixi run --environment py311 build-pandas --editable` |
| Run tests | Every code change, scoped to related tests (existing or added) | [contributing_codebase.rst](doc/source/development/contributing_codebase.rst) | `conda run -n pandas-dev pytest pandas/tests/<area>/test_<file>.py::test_name` | `pixi run --environment py311 python -m pytest pandas/tests/<area>/test_<file>.py::test_name` |
| Run doctests | A docstring `Examples` section changed | [contributing_docstring.rst](doc/source/development/contributing_docstring.rst) | `conda run -n pandas-dev env PANDAS_FUTURE_PYTHON_SCALARS=1 python -m pytest --doctest-modules --doctest-cython pandas/<module>.py` | `pixi run --environment doctests ci-doctests` |
| Linting / code style | Every change, on the changed files, before handing work back | [contributing_codebase.rst](doc/source/development/contributing_codebase.rst) | `conda run -n pandas-dev pre-commit run --files <changed files>` | `pixi run --environment typing pre-commit run --files <changed files>` |
| Validate typing | Type annotations or `.pyi` stubs changed | [contributing_codebase.rst](doc/source/development/contributing_codebase.rst) | `conda run -n pandas-dev pre-commit run --hook-stage manual mypy --files <changed files>` (also `pyright`, `stubtest`) | `pixi run --environment typing ci-typing` |
| Build documentation | Files under `doc/` changed | [contributing_documentation.rst](doc/source/development/contributing_documentation.rst) | `conda run -n pandas-dev python doc/make.py html` | `pixi run --environment docs-and-web ci-build-documentation html` |
| Run benchmarks | A change is performance-motivated | [contributing_codebase.rst](doc/source/development/contributing_codebase.rst) | `cd asv_bench && conda run -n pandas-dev asv continuous -f 1.1 upstream/main HEAD -b ^<pattern>` | `cd asv_bench && pixi run --environment asv asv continuous -f 1.1 upstream/main HEAD -b ^<pattern>` |

After the initial build, confirm a working local build with:

```bash
# conda
conda run -n pandas-dev python -c "import pandas as pd; print(pd.__version__)"
# pixi
pixi run --environment py311 python -c "import pandas as pd; print(pd.__version__)"
```

The version should resemble a local development build like `<version>.dev0+<integer>.<git hash>`, for example `3.1.0.dev0+880.g2b9e661fbb`.
If the import fails or the version looks like a released wheel, rerun a command in the "Build / rebuild" row.

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
- [contributing_docstring.rst](doc/source/development/contributing_docstring.rst) — docstring conventions
- [policies.rst](doc/source/development/policies.rst) — versioning and deprecation policy
