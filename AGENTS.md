# pandas Agent Instructions

## Project Overview
`pandas` is an open source, BSD-licensed library providing high-performance, easy-to-use data structures and data analysis tools for the Python programming language.

## Purpose
- Assist contributors by suggesting code changes, tests, and documentation edits for the pandas repository while preserving stability and compatibility.

## Persona & Tone
- Concise, neutral, code-focused. Prioritize correctness, readability, and tests.

## Project Guidelines
- Be sure to follow all guidelines for contributing to the codebase specified at https://pandas.pydata.org/docs/development/contributing_codebase.html
- These guidelines are also available in the following local files, which should be loaded into context and adhered to
    - doc/source/development/contributing_codebase.rst
    - doc/source/development/contributing_docstring.rst
    - doc/source/development/contributing_documentation.rst
    - doc/source/development/contributing.rst

## Decision heuristics
- Favor small, backward-compatible changes with tests.
- If a change would be breaking, propose it behind a deprecation path and document the rationale.
- Prefer readability over micro-optimizations unless benchmarks are requested.
- Add tests for behavioral changes; update docs only after code change is final.

## Type hints guidance (summary)
- Prefer PEP 484 style and types in pandas._typing when appropriate.
- Avoid unnecessary use of typing.cast; prefer refactors that convey types to type-checkers.
- Use builtin generics (list, dict) when possible.

## Docstring guidance (summary)
- Follow NumPy / numpydoc conventions used across the repo: short summary, extended summary, Parameters, Returns/Yields, See Also, Notes, Examples.
- Ensure examples are deterministic, import numpy/pandas as documented, and pass doctest rules used by docs validation.
- Preserve formatting rules: triple double-quotes, no blank line before/after docstring, parameter formatting ("name : type, default ..."), types and examples conventions.

## Pull Requests (summary)
- Pull request titles should be descriptive and include one of the following prefixes:
    - ENH: Enhancement, new functionality
    - BUG: Bug fix
    - DOC: Additions/updates to documentation
    - TST: Additions/updates to tests
    - BLD: Updates to the build process/scripts
    - PERF: Performance improvement
    - TYP: Type annotations
    - CLN: Code cleanup
- Pull request descriptions should follow the template, and **succinctly** describe the change being made. Usually a few sentences is sufficient.
- Pull requests which are resolving an existing Github Issue should include a link to the issue in the PR Description.
- Do not add summaries or additional comments to individual commit messages. The single PR description is sufficient.
