.. _contributing_code_guidelines:

===========================================
Pandas Code Style and Contribution Guide
===========================================

This single, authoritative, Sphinx-ready reference consolidates the project's
coding conventions, testing requirements, documentation practices, and
contribution workflow. It is designed to be concise, reviewer-friendly, and
actionable. Keep PRs small, provide tests, run pre-commit locally, and link
the related issue or discussion (see :issue: `33851`).

.. contents::
   :local:
   :depth: 2

Quick Start
===========

- Create a focused branch from upstream/main:
  ``git checkout -b doc/code-guidelines/phase-1-migrate-rst``.
- Add this file as ``doc/source/development/code_guidelines.rst`` and a one-line
  whatsnew note in ``doc/source/whatsnew/next.rst``.
- Run pre-commit hooks locally and add minimal, focused tests for code changes.
- Open a PR referencing issue #33851 and include the reviewer checklist below.

Philosophy
==========

- Clarity first: code should read like well-written prose; names and signatures
  must be intention-revealing.
- Safety and stability: avoid breaking user-facing APIs; prefer deprecation with
  clear migration guidance when necessary.
- Testability: every behavioral change must include tests demonstrating correct
  behavior for typical and edge cases.
- Reviewer empathy: small, well-documented changes with tests and a whatsnew
  entry are reviewed and merged faster.

Style and Formatting
====================

Line length
-----------
Maximum 88 characters. Use ``black`` to enforce formatting.

Imports
-------
Group imports in this order and alphabetize within groups:

1. Standard library
2. Third-party (e.g., ``numpy``, ``scipy``)
3. pandas
4. Local imports

Use ``isort`` via pre-commit to maintain ordering.

Naming
------
- Functions and variables: ``snake_case``
- Classes: ``CapWords``
- Constants: ``ALL_CAPS``

Docstrings
----------
Use the NumPy/SciPy (numpydoc) style for public APIs. Include a short summary,
parameters, returns, raises (if applicable), and concise examples when helpful.
Use Sphinx roles (``:func:``, ``:class:``, ``:meth:``) for cross-references.

Formatting tools
----------------
- ``black``: formatting
- ``ruff``: linting
- ``isort``: import sorting
- ``clang-format``: C/C++ style where applicable
- ``pre-commit``: orchestrates hooks

Idiomatic Python Guidelines
===========================

Context managers
----------------
Prefer ``with`` for resource management.

.. code-block:: python

    # Recommended
    with open("path") as f:
        data = f.read()

Boolean checks
--------------
Be explicit: use ``is None`` for None checks and avoid ambiguous truthiness.

.. code-block:: python

    # Recommended
    if value is None:
        ...

    if len(collection) == 0:
        ...

Exception handling
------------------
Catch specific exceptions only. Avoid bare ``except:`` or catching broad
``Exception`` without re-raising.

.. code-block:: python

    # Recommended
    try:
        x = int(s)
    except ValueError:
        x = float("nan")

Mutability and side effects
---------------------------
Prefer non-mutating APIs. If an API mutates an argument, document the behavior
explicitly and include tests asserting the mutation.

Type Hints and Annotations
==========================

- New public functions and significant refactors should include PEP 484 type
  annotations.
- Prefer pandas-specific aliases from ``pandas._typing`` (for example ``Dtype``).
- Avoid overusing ``cast``; refactor code to be type-checker-friendly where
  practical.
- Validate with ``mypy``/``pyright`` via pre-commit hooks.

.. code-block:: bash

    pre-commit run --hook-stage manual mypy --all-files
    pre-commit run --hook-stage manual pyright --all-files

Code Standards Checklist
========================

Before opening a PR, ensure:

- No trailing whitespace and no hard tabs.
- Code formatted with ``black``; imports sorted with ``isort``.
- No bare ``except`` clauses; catch specific exceptions only.
- Prefer Pythonic iteration over C-style index loops.
- Avoid unnecessary mutation of inputs.
- Tests added for new behavior; update existing tests if behavior changes.
- Docstrings updated for public APIs.
- Run pre-commit hooks locally.

Pre-commit and Local Validation
===============================

Install and run pre-commit locally to catch issues early.

.. code-block:: bash

    pip install pre-commit
    pre-commit install
    pre-commit run --files <modified-files>

To run all hooks between upstream/main and HEAD:

.. code-block:: bash

    pre-commit run --from-ref=upstream/main --to-ref=HEAD --all-files

Slow checks can be run manually:

.. code-block:: bash

    pre-commit run --hook-stage manual --all-files

Testing and Test-Driven Development
===================================

Testing philosophy
------------------
Adopt a tests-first mindset when feasible: write the minimal failing test,
implement the change, and iterate. Keep tests small, deterministic, and readable.

Test placement
--------------
Place tests in the relevant ``pandas/tests/...`` subdirectory. Prefer the most
specific location for the functionality being tested.

Common test locations (quick reference)
- tslibs helpers: ``tests.tslibs``, ``tests.scalar``
- C-extension libs: ``tests.libs``, ``tests.groupby.test_libgroupby``
- Arithmetic/comparison: ``tests.arithmetic``, ``tests.frame.test_arithmetic``
- Reductions: ``tests.reductions``, ``tests.frame.test_reductions``
- Indexing: ``tests.indexing``, ``tests.frame.test_getitem``, ``tests.series.test_getitem``
- IO: ``tests.io``
- Plotting: ``tests.plotting``
- ExtensionArray: ``tests.arrays``, shared EA tests in ``tests.extension``

Pytest idioms
-------------
- Use functional tests: ``def test_*``.
- Use bare ``assert`` for scalars.
- Use ``@pytest.mark.parametrize`` for multiple cases.
- Use ``pandas._testing`` utilities (``tm.assert_series_equal``, ``tm.assert_frame_equal``).
- Use ``tm.assert_produces_warning`` and ``pytest.raises(..., match="...")`` for warnings and exceptions.

Doctests and examples
---------------------
Ensure examples in docstrings are stable and deterministic. Use ``# doctest: +SKIP``
only when necessary for unstable examples.

Running tests locally
---------------------

.. code-block:: bash

    # Run a focused file
    pytest pandas/path/to/test_module.py -k "regex"

    # Parallelize locally (use with care)
    pytest -n auto pandas -m "not slow and not network and not db and not single_cpu"

    # Reproduce flaky behavior
    export PYTHONHASHSEED=314159265
    pytest ...

Property-based testing
----------------------
Use Hypothesis selectively for complex input domains. Prefer parametrized tests
for simple, enumerated cases to keep the suite fast.

Performance and Benchmarks
==========================

- Use ASV benchmarks in ``pandas/asv_bench`` for performance-sensitive changes.
- Target specific benchmarks with ``-b`` to keep runs tractable.
- Include a short benchmark summary in PRs when relevant.

Documentation and whatsnew
==========================

- Add a whatsnew entry under ``doc/source/whatsnew/vx.y.z.rst`` or ``next.rst``.
- Use Sphinx roles (``:func:``, ``:class:``, ``:meth:``) and the issue role for references.
- Keep whatsnew entries concise and include the issue/PR number.

Example whatsnew entry
----------------------

.. code-block:: rst

    Bug fixes
    ---------
    * Fixed ordering bug in :meth:`pandas.DataFrame.foo` (:issue:`12345`).

Pull Request Best Practices
===========================

Crafting the PR
---------------
- Use a title prefix: ENH, BUG, DOC, TST, PERF, TYP, CLN.
- Provide a 2–4 sentence summary, link to the issue, and list tests added.
- Include the reviewer checklist in the PR body.

Reviewer checklist (include in PR body)
- [ ] Summary and motivation provided
- [ ] Link to issue
- [ ] Tests added and passing locally
- [ ] Docs/whatsnew updated
- [ ] Pre-commit hooks run locally
- [ ] Benchmark evidence included (if applicable)

Commit hygiene
--------------
- One logical change per commit.
- Use imperative, present-tense commit messages (e.g., "Fix X", "Add Y tests").

CI Checks, Common Failures, and Remediation
===========================================

Common CI failures
- Formatting / import ordering: run ``black`` and ``isort`` locally.
- Lint failures (``ruff``): fix unused imports and simple style issues.
- Type-check failures: add concise annotations or refactor code.
- Docstring/doctest failures: make examples deterministic or mark skips.
- Test failures: run the focused tests locally and stabilize flaky behavior.

Typical remediation commands
----------------------------

.. code-block:: bash

    pre-commit run --hook-stage manual black isort ruff --all-files
    pre-commit run --hook-stage manual mypy --all-files
    make -C doc html
    pytest pandas/path/to/test_module.py -k "your_test_regex"

Debugging CI failures
---------------------
- Reproduce CI failures locally where possible (use CI's environment info).
- Isolate failing tests with ``-k`` and run single-threaded.
- Match Python and dependency versions when reproducing type or runtime errors.

Contribution Workflow and Etiquette
===================================

Before you start
- Search open issues to avoid duplicated effort.
- Claim an unassigned issue with a short plan and ETA.
- For large changes, open a design discussion (issue or draft PR).

Branching and syncing
- Keep local ``main`` synced: ``git fetch upstream && git checkout main && git merge upstream/main --ff-only``.
- Create focused branches: ``git checkout -b feature/concise-title``.

Responding to review
- Address comments promptly and in small commits.
- If you disagree, respond with rationale and propose alternatives.
- Update tests and whatsnew when behavior changes.

Security and sensitive data
- Do not commit secrets. If leaked, rotate credentials immediately and notify maintainers.
- Report security issues via the repository's security disclosure process.

Repository Layout (developer view)
==================================
- ``doc/``: documentation sources
- ``doc/source/development/``: developer-facing docs (this file)
- ``pandas/``: library source
- ``pandas/tests/``: test suite
- ``ci/``: CI scripts and validators
- ``pandas/asv_bench/``: performance benchmarks
- ``.pre-commit-config.yaml``: hook configuration
- ``pyproject.toml``: linting/format defaults

Phase 2: CI and Linter Rule Extraction (Roadmap)
================================================

Goal
----
Translate CI and pre-commit checks into concise, human-readable rules with:

- A one-line requirement
- Minimal bad/good examples
- A "Derived from" pointer (script and approximate line range)
- A short remediation note

Approach
--------
1. Parse ``.pre-commit-config.yaml``, ``pyproject.toml``, and ``ci/`` scripts.
2. Extract user-facing messages and map them to rules.
3. Prioritize high-frequency CI failures (formatting, imports, type checks).
4. Author small PRs (e.g., top 10 rules) with examples and derived-from citations.

Example rule entry (template)
-----------------------------

.. code-block:: rst

    .. _rule-isort:

    Import ordering (isort)
    -----------------------
    Requirement
      Group imports: standard library; third-party; pandas; local. Alphabetize within groups.
    Bad example
      ``import pandas as pd\nimport numpy as np\nfrom os import path``
    Good example
      ``from os import path\nimport numpy as np\nimport pandas as pd``
    How to check locally
      ``pre-commit run isort --files <modified-files>``
    Derived from
      ``.pre-commit-config.yaml`` (isort hook)

Onboarding Checklist (first-time contributors)
==============================================

1. Fork the repository on GitHub.
2. Sparse-clone if bandwidth is limited (see Appendix).
3. Create a feature branch from upstream/main.
4. Copy this file into ``doc/source/development/``.
5. Add a one-line whatsnew entry in ``doc/source/whatsnew/next.rst``.
6. Run pre-commit locally and add tests.
7. Commit, push to your fork, and open a PR referencing :issue: `33851`.

Appendix: Commands and Recipes
==============================

Sparse clone (low bandwidth)

.. code-block:: bash

    git clone --depth=1 --filter=blob:none --no-checkout https://github.com/pandas-dev/pandas.git pandas-sparse
    cd pandas-sparse
    git sparse-checkout init --cone
    git sparse-checkout set doc doc/source/whatsnew ci .pre-commit-config.yaml pyproject.toml
    git checkout main

Branch and files

.. code-block:: bash

    git checkout -b doc/code-guidelines/phase-1-migrate-rst
    mkdir -p doc/source/development
    # paste this file into doc/source/development/code_guidelines.rst
    # add one-line whatsnew into doc/source/whatsnew/next.rst
    git add doc/source/development/code_guidelines.rst doc/source/whatsnew/next.rst
    git commit -m "DOC: add consolidated code guidelines (phase 1) #33851"

Create a patch (if you cannot push)

.. code-block:: bash

    git format-patch -1 HEAD --stdout > phase-1-migrate-rst.patch
    # Maintainers can apply with:
    # git am < phase-1-migrate-rst.patch

Resume and visibility guidance
==============================

- Keep PRs small and document impact clearly. A merged, well-documented PR is a
  high-signal contribution.
- After merge, add a concise resume bullet with the PR link and a short outcome.

Acknowledgements
================
This file consolidates developer guidance into a durable reference. Phase 2 will
add CI-derived exact rules and minimal bad/good examples with script citations
to complete the canonical source of truth.
