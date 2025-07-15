# This file is part of Hypothesis, which may be found at
# https://github.com/HypothesisWorks/hypothesis/
#
# Copyright the Hypothesis Authors.
# Individual contributors are listed in AUTHORS.rst and the git log.
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at https://mozilla.org/MPL/2.0/.

"""The settings module configures runtime options for Hypothesis.

Either an explicit settings object can be used or the default object on
this module can be modified.
"""

import contextlib
import datetime
import inspect
import os
import warnings
from collections.abc import Collection, Generator, Sequence
from enum import Enum, EnumMeta, IntEnum, unique
from typing import (
    TYPE_CHECKING,
    Any,
    ClassVar,
    Optional,
    TypeVar,
    Union,
)

from hypothesis.errors import (
    HypothesisDeprecationWarning,
    InvalidArgument,
)
from hypothesis.internal.conjecture.providers import AVAILABLE_PROVIDERS
from hypothesis.internal.reflection import get_pretty_function_description
from hypothesis.internal.validation import check_type, try_convert
from hypothesis.utils.conventions import not_set
from hypothesis.utils.dynamicvariables import DynamicVariable

if TYPE_CHECKING:
    from hypothesis.database import ExampleDatabase

__all__ = ["settings"]

T = TypeVar("T")
all_settings: list[str] = [
    "max_examples",
    "derandomize",
    "database",
    "verbosity",
    "phases",
    "stateful_step_count",
    "report_multiple_bugs",
    "suppress_health_check",
    "deadline",
    "print_blob",
    "backend",
]


@unique
class Verbosity(IntEnum):
    """Options for the |settings.verbosity| argument to |@settings|."""

    quiet = 0
    """
    Hypothesis will not print any output, not even the final falsifying example.
    """

    normal = 1
    """
    Standard verbosity. Hypothesis will print the falsifying example, alongside
    any notes made with |note| (only for the falsfying example).
    """

    verbose = 2
    """
    Increased verbosity. In addition to everything in |Verbosity.normal|, Hypothesis
    will print each example as it tries it, as well as any notes made with |note|
    for every example. Hypothesis will also print shrinking attempts.
    """

    debug = 3
    """
    Even more verbosity. Useful for debugging Hypothesis internals. You probably
    don't want this.
    """

    def __repr__(self) -> str:
        return f"Verbosity.{self.name}"


@unique
class Phase(IntEnum):
    """Options for the |settings.phases| argument to |@settings|."""

    explicit = 0
    """
    Controls whether explicit examples are run.
    """

    reuse = 1
    """
    Controls whether previous examples will be reused.
    """

    generate = 2
    """
    Controls whether new examples will be generated.
    """

    target = 3
    """
    Controls whether examples will be mutated for targeting.
    """

    shrink = 4
    """
    Controls whether examples will be shrunk.
    """

    explain = 5
    """
    Controls whether Hypothesis attempts to explain test failures.

    The explain phase has two parts, each of which is best-effort - if Hypothesis
    can't find a useful explanation, we'll just print the minimal failing example.
    """

    def __repr__(self) -> str:
        return f"Phase.{self.name}"


class HealthCheckMeta(EnumMeta):
    def __iter__(self):
        deprecated = (HealthCheck.return_value, HealthCheck.not_a_test_method)
        return iter(x for x in super().__iter__() if x not in deprecated)


@unique
class HealthCheck(Enum, metaclass=HealthCheckMeta):
    """Arguments for :attr:`~hypothesis.settings.suppress_health_check`.

    Each member of this enum is a specific health check to suppress.

    Hypothesis' health checks are designed to detect and warn you about performance
    problems where your tests are slow, inefficient, or generating very large examples.

    If this is expected, e.g. when generating large arrays or dataframes, you can selectively
    disable them with the :obj:`~hypothesis.settings.suppress_health_check` setting.
    The argument for this parameter is a list with elements drawn from any of
    the class-level attributes of the HealthCheck class.
    Using a value of ``list(HealthCheck)`` will disable all health checks.
    """

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}.{self.name}"

    @classmethod
    def all(cls) -> list["HealthCheck"]:
        # Skipping of deprecated attributes is handled in HealthCheckMeta.__iter__
        note_deprecation(
            "`HealthCheck.all()` is deprecated; use `list(HealthCheck)` instead.",
            since="2023-04-16",
            has_codemod=True,
            stacklevel=1,
        )
        return list(HealthCheck)

    data_too_large = 1
    """Checks if too many examples are aborted for being too large.

    This is measured by the number of random choices that Hypothesis makes
    in order to generate something, not the size of the generated object.
    For example, choosing a 100MB object from a predefined list would take
    only a few bits, while generating 10KB of JSON from scratch might trigger
    this health check.
    """

    filter_too_much = 2
    """Check for when the test is filtering out too many examples, either
    through use of |assume| or |.filter|, or occasionally for Hypothesis
    internal reasons."""

    too_slow = 3
    """Check for when your data generation is extremely slow and likely to hurt
    testing."""

    return_value = 5
    """Deprecated; we always error if a test returns a non-None value."""

    large_base_example = 7
    """Checks if the natural example to shrink towards is very large."""

    not_a_test_method = 8
    """Deprecated; we always error if |@given| is applied
    to a method defined by :class:`python:unittest.TestCase` (i.e. not a test)."""

    function_scoped_fixture = 9
    """Checks if |@given| has been applied to a test
    with a pytest function-scoped fixture. Function-scoped fixtures run once
    for the whole function, not once per example, and this is usually not what
    you want.

    Because of this limitation, tests that need to set up or reset
    state for every example need to do so manually within the test itself,
    typically using an appropriate context manager.

    Suppress this health check only in the rare case that you are using a
    function-scoped fixture that does not need to be reset between individual
    examples, but for some reason you cannot use a wider fixture scope
    (e.g. session scope, module scope, class scope).

    This check requires the :ref:`Hypothesis pytest plugin<pytest-plugin>`,
    which is enabled by default when running Hypothesis inside pytest."""

    differing_executors = 10
    """Checks if |@given| has been applied to a test
    which is executed by different :ref:`executors<custom-function-execution>`.
    If your test function is defined as a method on a class, that class will be
    your executor, and subclasses executing an inherited test is a common way
    for things to go wrong.

    The correct fix is often to bring the executor instance under the control
    of hypothesis by explicit parametrization over, or sampling from,
    subclasses, or to refactor so that |@given| is
    specified on leaf subclasses."""

    nested_given = 11
    """Checks if |@given| is used inside another
    |@given|. This results in quadratic generation and
    shrinking behavior, and can usually be expressed more cleanly by using
    :func:`~hypothesis.strategies.data` to replace the inner
    |@given|.

    Nesting @given can be appropriate if you set appropriate limits for the
    quadratic behavior and cannot easily reexpress the inner function with
    :func:`~hypothesis.strategies.data`. To suppress this health check, set
    ``suppress_health_check=[HealthCheck.nested_given]`` on the outer
    |@given|. Setting it on the inner
    |@given| has no effect. If you have more than one
    level of nesting, add a suppression for this health check to every
    |@given| except the innermost one.
    """


class duration(datetime.timedelta):
    """A timedelta specifically measured in milliseconds."""

    def __repr__(self) -> str:
        ms = self.total_seconds() * 1000
        return f"timedelta(milliseconds={int(ms) if ms == int(ms) else ms!r})"


def is_in_ci() -> bool:
    # GitHub Actions, Travis CI and AppVeyor have "CI"
    # Azure Pipelines has "TF_BUILD"
    # GitLab CI has "GITLAB_CI"
    return "CI" in os.environ or "TF_BUILD" in os.environ or "GITLAB_CI" in os.environ


default_variable = DynamicVariable[Optional["settings"]](None)


def _validate_choices(name: str, value: T, *, choices: Sequence[object]) -> T:
    if value not in choices:
        msg = f"Invalid {name}, {value!r}. Valid choices: {choices!r}"
        raise InvalidArgument(msg)
    return value


def _validate_max_examples(max_examples: int) -> int:
    check_type(int, max_examples, name="max_examples")
    if max_examples < 1:
        raise InvalidArgument(
            f"max_examples={max_examples!r} must be at least one. If you want "
            "to disable generation entirely, use phases=[Phase.explicit] instead."
        )
    return max_examples


def _validate_database(
    database: Optional["ExampleDatabase"],
) -> Optional["ExampleDatabase"]:
    from hypothesis.database import ExampleDatabase

    if database is None or isinstance(database, ExampleDatabase):
        return database
    raise InvalidArgument(
        "Arguments to the database setting must be None or an instance of "
        "ExampleDatabase. Use one of the database classes in "
        "hypothesis.database"
    )


def _validate_phases(phases: Collection[Phase]) -> Sequence[Phase]:
    phases = tuple(phases)
    for phase in phases:
        if not isinstance(phase, Phase):
            raise InvalidArgument(f"{phase!r} is not a valid phase")
    return tuple(phase for phase in list(Phase) if phase in phases)


def _validate_stateful_step_count(stateful_step_count: int) -> int:
    check_type(int, stateful_step_count, name="stateful_step_count")
    if stateful_step_count < 1:
        raise InvalidArgument(
            f"stateful_step_count={stateful_step_count!r} must be at least one."
        )
    return stateful_step_count


def _validate_suppress_health_check(suppressions):
    suppressions = try_convert(tuple, suppressions, "suppress_health_check")
    for health_check in suppressions:
        if not isinstance(health_check, HealthCheck):
            raise InvalidArgument(
                f"Non-HealthCheck value {health_check!r} of type {type(health_check).__name__} "
                "is invalid in suppress_health_check."
            )
        if health_check in (HealthCheck.return_value, HealthCheck.not_a_test_method):
            note_deprecation(
                f"The {health_check.name} health check is deprecated, because this is always an error.",
                since="2023-03-15",
                has_codemod=False,
                stacklevel=2,
            )
    return suppressions


def _validate_deadline(
    x: Union[int, float, datetime.timedelta, None],
) -> Optional[duration]:
    if x is None:
        return x
    invalid_deadline_error = InvalidArgument(
        f"deadline={x!r} (type {type(x).__name__}) must be a timedelta object, "
        "an integer or float number of milliseconds, or None to disable the "
        "per-test-case deadline."
    )
    if isinstance(x, (int, float)):
        if isinstance(x, bool):
            raise invalid_deadline_error
        try:
            x = duration(milliseconds=x)
        except OverflowError:
            raise InvalidArgument(
                f"deadline={x!r} is invalid, because it is too large to represent "
                "as a timedelta. Use deadline=None to disable deadlines."
            ) from None
    if isinstance(x, datetime.timedelta):
        if x <= datetime.timedelta(0):
            raise InvalidArgument(
                f"deadline={x!r} is invalid, because it is impossible to meet a "
                "deadline <= 0. Use deadline=None to disable deadlines."
            )
        return duration(seconds=x.total_seconds())
    raise invalid_deadline_error


def _validate_backend(backend: str) -> str:
    if backend not in AVAILABLE_PROVIDERS:
        if backend == "crosshair":  # pragma: no cover
            install = '`pip install "hypothesis[crosshair]"` and try again.'
            raise InvalidArgument(f"backend={backend!r} is not available.  {install}")
        raise InvalidArgument(
            f"backend={backend!r} is not available - maybe you need to install a plugin?"
            f"\n    Installed backends: {sorted(AVAILABLE_PROVIDERS)!r}"
        )
    return backend


class settingsMeta(type):
    def __init__(cls, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @property
    def default(cls) -> Optional["settings"]:
        v = default_variable.value
        if v is not None:
            return v
        if getattr(settings, "_current_profile", None) is not None:
            assert settings._current_profile is not None
            settings.load_profile(settings._current_profile)
            assert default_variable.value is not None
        return default_variable.value

    def __setattr__(cls, name: str, value: object) -> None:
        if name == "default":
            raise AttributeError(
                "Cannot assign to the property settings.default - "
                "consider using settings.load_profile instead."
            )
        elif not name.startswith("_"):
            raise AttributeError(
                f"Cannot assign hypothesis.settings.{name}={value!r} - the settings "
                "class is immutable.  You can change the global default "
                "settings with settings.load_profile, or use @settings(...) "
                "to decorate your test instead."
            )
        super().__setattr__(name, value)

    def __repr__(cls):
        return "hypothesis.settings"


class settings(metaclass=settingsMeta):
    """
    A settings object controls the following aspects of test behavior:
    |~settings.max_examples|, |~settings.derandomize|, |~settings.database|,
    |~settings.verbosity|, |~settings.phases|, |~settings.stateful_step_count|,
    |~settings.report_multiple_bugs|, |~settings.suppress_health_check|,
    |~settings.deadline|, |~settings.print_blob|, and |~settings.backend|.

    A settings object can be applied as a decorator to a test function, in which
    case that test function will use those settings. A test may only have one
    settings object applied to it. A settings object can also be passed to
    |settings.register_profile| or as a parent to another |settings|.

    Attribute inheritance
    ---------------------

    Settings objects are immutable once created. When a settings object is created,
    it uses the value specified for each attribute. Any attribute which is
    not specified will inherit from its value in the ``parent`` settings object.
    If ``parent`` is not passed, any attributes which are not specified will inherit
    from the currently active settings profile instead.

    For instance, ``settings(max_examples=10)`` will have a ``max_examples`` of ``10``,
    and the value of all other attributes will be equal to its value in the
    currently active settings profile.

    A settings object is immutable once created. Changes made from activating a new
    settings profile with |settings.load_profile| will be reflected in
    settings objects created after the profile was made active, but not in existing
    settings objects.

    Built-in profiles
    -----------------

    While you can register additional profiles with |settings.register_profile|,
    Hypothesis comes with two built-in profiles: ``default`` and ``ci``.

    The ``default`` profile is active by default, unless one of the ``CI``,
    ``TF_BUILD``, or ``GITLAB_CI`` environment variables are set (to any value),
    in which case the ``CI`` profile will be active by default.

    The attributes of the currently active settings profile can be retrieved with
    ``settings()`` (so ``settings().max_examples`` is the currently active default
    for |settings.max_examples|).

    The settings attributes for the built-in profiles are as follows:

    .. code-block:: python

        default = settings.register_profile(
            "default",
            max_examples=100,
            derandomize=False,
            database=not_set,  # see settings.database for details
            verbosity=Verbosity.normal,
            phases=tuple(Phase),
            stateful_step_count=50,
            report_multiple_bugs=True,
            suppress_health_check=(),
            deadline=duration(milliseconds=200),
            print_blob=False,
            backend="hypothesis",
        )

        ci = settings.register_profile(
            "ci",
            parent=default,
            derandomize=True,
            deadline=None,
            database=None,
            print_blob=True,
            suppress_health_check=[HealthCheck.too_slow],
        )

    You can configure either of the built-in profiles with |settings.register_profile|:

    .. code-block:: python

        # run more examples in CI
        settings.register_profile(
            "ci",
            settings.get_profile("ci"),
            max_examples=1000,
        )
    """

    _profiles: ClassVar[dict[str, "settings"]] = {}
    _current_profile: ClassVar[Optional[str]] = None

    def __init__(
        self,
        parent: Optional["settings"] = None,
        *,
        # This looks pretty strange, but there's good reason: we want Mypy to detect
        # bad calls downstream, but not to freak out about the `= not_set` part even
        # though it's not semantically valid to pass that as an argument value.
        # The intended use is "like **kwargs, but more tractable for tooling".
        max_examples: int = not_set,  # type: ignore
        derandomize: bool = not_set,  # type: ignore
        database: Optional["ExampleDatabase"] = not_set,  # type: ignore
        verbosity: "Verbosity" = not_set,  # type: ignore
        phases: Collection["Phase"] = not_set,  # type: ignore
        stateful_step_count: int = not_set,  # type: ignore
        report_multiple_bugs: bool = not_set,  # type: ignore
        suppress_health_check: Collection["HealthCheck"] = not_set,  # type: ignore
        deadline: Union[int, float, datetime.timedelta, None] = not_set,  # type: ignore
        print_blob: bool = not_set,  # type: ignore
        backend: str = not_set,  # type: ignore
    ) -> None:
        self._in_definition = True

        if parent is not None:
            check_type(settings, parent, "parent")
        if derandomize not in (not_set, False):
            if database not in (not_set, None):  # type: ignore
                raise InvalidArgument(
                    "derandomize=True implies database=None, so passing "
                    f"{database=} too is invalid."
                )
            database = None

        # fallback is None if we're creating the default settings object, and
        # the parent (or default settings object) otherwise
        self._fallback = parent or settings.default
        self._max_examples = (
            self._fallback.max_examples  # type: ignore
            if max_examples is not_set  # type: ignore
            else _validate_max_examples(max_examples)
        )
        self._derandomize = (
            self._fallback.derandomize  # type: ignore
            if derandomize is not_set  # type: ignore
            else _validate_choices("derandomize", derandomize, choices=[True, False])
        )
        if database is not not_set:  # type: ignore
            database = _validate_database(database)
        self._database = database
        self._cached_database = None
        self._verbosity = (
            self._fallback.verbosity  # type: ignore
            if verbosity is not_set  # type: ignore
            else _validate_choices("verbosity", verbosity, choices=tuple(Verbosity))
        )
        self._phases = (
            self._fallback.phases  # type: ignore
            if phases is not_set  # type: ignore
            else _validate_phases(phases)
        )
        self._stateful_step_count = (
            self._fallback.stateful_step_count  # type: ignore
            if stateful_step_count is not_set  # type: ignore
            else _validate_stateful_step_count(stateful_step_count)
        )
        self._report_multiple_bugs = (
            self._fallback.report_multiple_bugs  # type: ignore
            if report_multiple_bugs is not_set  # type: ignore
            else _validate_choices(
                "report_multiple_bugs", report_multiple_bugs, choices=[True, False]
            )
        )
        self._suppress_health_check = (
            self._fallback.suppress_health_check  # type: ignore
            if suppress_health_check is not_set  # type: ignore
            else _validate_suppress_health_check(suppress_health_check)
        )
        self._deadline = (
            self._fallback.deadline  # type: ignore
            if deadline is not_set
            else _validate_deadline(deadline)
        )
        self._print_blob = (
            self._fallback.print_blob  # type: ignore
            if print_blob is not_set  # type: ignore
            else _validate_choices("print_blob", print_blob, choices=[True, False])
        )
        self._backend = (
            self._fallback.backend  # type: ignore
            if backend is not_set  # type: ignore
            else _validate_backend(backend)
        )

        self._in_definition = False

    @property
    def max_examples(self):
        """
        Once this many satisfying examples have been considered without finding any
        counter-example, Hypothesis will stop looking.

        Note that we might call your test function fewer times if we find a bug early
        or can tell that we've exhausted the search space; or more if we discard some
        examples due to use of .filter(), assume(), or a few other things that can
        prevent the test case from completing successfully.

        The default value is chosen to suit a workflow where the test will be part of
        a suite that is regularly executed locally or on a CI server, balancing total
        running time against the chance of missing a bug.

        If you are writing one-off tests, running tens of thousands of examples is
        quite reasonable as Hypothesis may miss uncommon bugs with default settings.
        For very complex code, we have observed Hypothesis finding novel bugs after
        *several million* examples while testing :pypi:`SymPy <sympy>`.
        If you are running more than 100k examples for a test, consider using our
        :ref:`integration for coverage-guided fuzzing <fuzz_one_input>` - it really
        shines when given minutes or hours to run.

        The default max examples is ``100``.
        """
        return self._max_examples

    @property
    def derandomize(self):
        """
        If True, seed Hypothesis' random number generator using a hash of the test
        function, so that every run will test the same set of examples until you
        update Hypothesis, Python, or the test function.

        This allows you to `check for regressions and look for bugs
        <https://blog.nelhage.com/post/two-kinds-of-testing/>`__ using separate
        settings profiles - for example running
        quick deterministic tests on every commit, and a longer non-deterministic
        nightly testing run.

        The default is ``False``. If running on CI, the default is ``True`` instead.
        """
        return self._derandomize

    @property
    def database(self):
        """
        An instance of |ExampleDatabase| that will be used to save examples to
        and load previous examples from.

        If not set, a |DirectoryBasedExampleDatabase| is created in the current
        working directory under ``.hypothesis/examples``. If this location is
        unusable, e.g. due to the lack of read or write permissions, Hypothesis
        will emit a warning and fall back to an |InMemoryExampleDatabase|.

        If ``None``, no storage will be used.

        See the :ref:`database documentation <database>` for a list of database
        classes, and how to define custom database classes.
        """
        from hypothesis.database import _db_for_path

        # settings.database has two conflicting requirements:
        # * The default settings should respect changes to set_hypothesis_home_dir
        #   in-between accesses
        # * `s.database is s.database` should be true, except for the default settings
        #
        # We therefore cache s.database for everything except the default settings,
        # which always recomputes dynamically.
        if self._fallback is None:
            # if self._fallback is None, we are the default settings, at which point
            # we should recompute the database dynamically
            assert self._database is not_set
            return _db_for_path(not_set)

        # otherwise, we cache the database
        if self._cached_database is None:
            self._cached_database = (
                self._fallback.database if self._database is not_set else self._database
            )
        return self._cached_database

    @property
    def verbosity(self):
        """
        Control the verbosity level of Hypothesis messages.

        To see what's going on while Hypothesis runs your tests, you can turn
        up the verbosity setting.

        .. code-block:: pycon

            >>> from hypothesis import settings, Verbosity
            >>> from hypothesis.strategies import lists, integers
            >>> @given(lists(integers()))
            ... @settings(verbosity=Verbosity.verbose)
            ... def f(x):
            ...     assert not any(x)
            ... f()
            Trying example: []
            Falsifying example: [-1198601713, -67, 116, -29578]
            Shrunk example to [-1198601713]
            Shrunk example to [-128]
            Shrunk example to [32]
            Shrunk example to [1]
            [1]

        The four levels are |Verbosity.quiet|, |Verbosity.normal|,
        |Verbosity.verbose|, and |Verbosity.debug|. |Verbosity.normal| is the
        default. For |Verbosity.quiet|, Hypothesis will not print anything out,
        not even the final falsifying example. |Verbosity.debug| is basically
        |Verbosity.verbose| but a bit more so. You probably don't want it.

        If you are using :pypi:`pytest`, you may also need to :doc:`disable
        output capturing for passing tests <pytest:how-to/capture-stdout-stderr>`
        to see verbose output as tests run.
        """
        return self._verbosity

    @property
    def phases(self):
        """
        Control which phases should be run.

        Hypothesis divides tests into logically distinct phases.

        - |Phase.explicit|: Running explicit examples from |@example|.
        - |Phase.reuse|: Running examples from the database which previously failed.
        - |Phase.generate|: Generating new random examples.
        - |Phase.target|: Mutating examples for :ref:`targeted property-based
          testing <targeted>`. Requires |Phase.generate|.
        - |Phase.shrink|: Shrinking failing examples.
        - |Phase.explain|: Attempting to explain why a failure occurred.
          Requires |Phase.shrink|.

        Following the first failure, Hypothesis will (usually, depending on
        which |Phase| is enabled) track which lines of code are always run on
        failing but never on passing inputs. On 3.12+, this uses
        :mod:`sys.monitoring`, while 3.11 and earlier uses :func:`python:sys.settrace`.
        For python 3.11 and earlier, we therefore automatically disable the explain
        phase on PyPy, or if you are using :pypi:`coverage` or a debugger. If
        there are no clearly suspicious lines of code, :pep:`we refuse the
        temptation to guess <20>`.

        After shrinking to a minimal failing example, Hypothesis will try to find
        parts of the example -- e.g. separate args to |@given|
        -- which can vary freely without changing the result
        of that minimal failing example. If the automated experiments run without
        finding a passing variation, we leave a comment in the final report:

        .. code-block:: python

            test_x_divided_by_y(
                x=0,  # or any other generated value
                y=0,
            )

        Just remember that the *lack* of an explanation sometimes just means that
        Hypothesis couldn't efficiently find one, not that no explanation (or
        simpler failing example) exists.


        The phases setting provides you with fine grained control over which of
        these run, with each phase corresponding to a value on the |Phase| enum.

        The phases argument accepts a collection with any subset of these. e.g.
        ``settings(phases=[Phase.generate, Phase.shrink])`` will generate new examples
        and shrink them, but will not run explicit examples or reuse previous failures,
        while ``settings(phases=[Phase.explicit])`` will only run the explicit
        examples.
        """

        return self._phases

    @property
    def stateful_step_count(self):
        """
        The maximum number of times to call an additional |@rule| method in
        :ref:`stateful testing <stateful>` before we give up on finding a bug.

        Note that this setting is effectively multiplicative with max_examples,
        as each example will run for a maximum of ``stateful_step_count`` steps.

        The default stateful step count is ``50``.
        """
        return self._stateful_step_count

    @property
    def report_multiple_bugs(self):
        """
        Because Hypothesis runs the test many times, it can sometimes find multiple
        bugs in a single run.  Reporting all of them at once is usually very useful,
        but replacing the exceptions can occasionally clash with debuggers.
        If disabled, only the exception with the smallest minimal example is raised.

        The default value is ``True``.
        """
        return self._report_multiple_bugs

    @property
    def suppress_health_check(self):
        """
        A list of |HealthCheck| items to disable.
        """
        return self._suppress_health_check

    @property
    def deadline(self):
        """
        The maximum allowed duration of an individual test case, in milliseconds.
        You can pass an integer, float, or timedelta. If ``None``, the deadline
        is disabled entirely.

        We treat the deadline as a soft limit in some cases, where that would
        avoid flakiness due to timing variability.

        The default deadline is 200 milliseconds. If running on CI, the default is
        ``None`` instead.
        """
        return self._deadline

    @property
    def print_blob(self):
        """
        If set to ``True``, Hypothesis will print code for failing examples that
        can be used with |@reproduce_failure| to reproduce the failing example.

        The default value is ``False``. If running on CI, the default is ``True`` instead.
        """
        return self._print_blob

    @property
    def backend(self):
        """
        .. warning::

            EXPERIMENTAL AND UNSTABLE - see :ref:`alternative-backends`.

        The importable name of a backend which Hypothesis should use to generate
        primitive types. We support heuristic-random, solver-based, and fuzzing-based
        backends.
        """
        return self._backend

    def __call__(self, test: T) -> T:
        """Make the settings object (self) an attribute of the test.

        The settings are later discovered by looking them up on the test itself.
        """
        # Aliasing as Any avoids mypy errors (attr-defined) when accessing and
        # setting custom attributes on the decorated function or class.
        _test: Any = test

        # Using the alias here avoids a mypy error (return-value) later when
        # ``test`` is returned, because this check results in type refinement.
        if not callable(_test):
            raise InvalidArgument(
                "settings objects can be called as a decorator with @given, "
                f"but decorated {test=} is not callable."
            )
        if inspect.isclass(test):
            from hypothesis.stateful import RuleBasedStateMachine

            if issubclass(_test, RuleBasedStateMachine):
                attr_name = "_hypothesis_internal_settings_applied"
                if getattr(test, attr_name, False):
                    raise InvalidArgument(
                        "Applying the @settings decorator twice would "
                        "overwrite the first version; merge their arguments "
                        "instead."
                    )
                setattr(test, attr_name, True)
                _test.TestCase.settings = self
                return test  # type: ignore
            else:
                raise InvalidArgument(
                    "@settings(...) can only be used as a decorator on "
                    "functions, or on subclasses of RuleBasedStateMachine."
                )
        if hasattr(_test, "_hypothesis_internal_settings_applied"):
            # Can't use _hypothesis_internal_use_settings as an indicator that
            # @settings was applied, because @given also assigns that attribute.
            descr = get_pretty_function_description(test)
            raise InvalidArgument(
                f"{descr} has already been decorated with a settings object.\n"
                f"    Previous:  {_test._hypothesis_internal_use_settings!r}\n"
                f"    This:  {self!r}"
            )

        _test._hypothesis_internal_use_settings = self
        _test._hypothesis_internal_settings_applied = True
        return test

    def __setattr__(self, name: str, value: object) -> None:
        if not name.startswith("_") and not self._in_definition:
            raise AttributeError("settings objects are immutable")
        return super().__setattr__(name, value)

    def __repr__(self) -> str:
        bits = sorted(
            f"{name}={getattr(self, name)!r}"
            for name in all_settings
            if (name != "backend" or len(AVAILABLE_PROVIDERS) > 1)  # experimental
        )
        return "settings({})".format(", ".join(bits))

    def show_changed(self) -> str:
        bits = []
        for name in all_settings:
            value = getattr(self, name)
            if value != getattr(default, name):
                bits.append(f"{name}={value!r}")
        return ", ".join(sorted(bits, key=len))

    @staticmethod
    def register_profile(
        name: str,
        parent: Optional["settings"] = None,
        **kwargs: Any,
    ) -> None:
        """
        Register a settings object as a settings profile, under the name ``name``.
        The ``parent`` and ``kwargs`` arguments to this method are as for
        |settings|.

        If a settings profile already exists under ``name``, it will be overwritten.
        Registering a profile with the same name as the currently active profile
        will cause those changes to take effect in the active profile immediately,
        and do not require reloading the profile.

        Registered settings profiles can be retrieved later by name with
        |settings.get_profile|.
        """
        check_type(str, name, "name")
        # if we just pass the parent and no kwargs, like
        #   settings.register_profile(settings(max_examples=10))
        # then optimize out the pointless intermediate settings object which
        # would just forward everything to the parent.
        settings._profiles[name] = (
            parent
            if parent is not None and not kwargs
            else settings(parent=parent, **kwargs)
        )
        if settings._current_profile == name:
            settings.load_profile(name)

    @staticmethod
    def get_profile(name: str) -> "settings":
        """
        Returns the settings profile registered under ``name``. If no settings
        profile is registered under ``name``, raises |InvalidArgument|.
        """
        check_type(str, name, "name")
        try:
            return settings._profiles[name]
        except KeyError:
            raise InvalidArgument(f"Profile {name!r} is not registered") from None

    @staticmethod
    def load_profile(name: str) -> None:
        """
        Makes the settings profile registered under ``name`` the active profile.

        If no settings profile is registered under ``name``, raises |InvalidArgument|.
        """
        check_type(str, name, "name")
        settings._current_profile = name
        default_variable.value = settings.get_profile(name)


@contextlib.contextmanager
def local_settings(s: settings) -> Generator[settings, None, None]:
    with default_variable.with_value(s):
        yield s


def note_deprecation(
    message: str, *, since: str, has_codemod: bool, stacklevel: int = 0
) -> None:
    if since != "RELEASEDAY":
        date = datetime.date.fromisoformat(since)
        assert datetime.date(2021, 1, 1) <= date
    if has_codemod:
        message += (
            "\n    The `hypothesis codemod` command-line tool can automatically "
            "refactor your code to fix this warning."
        )
    warnings.warn(HypothesisDeprecationWarning(message), stacklevel=2 + stacklevel)


default = settings(
    max_examples=100,
    derandomize=False,
    database=not_set,  # type: ignore
    verbosity=Verbosity.normal,
    phases=tuple(Phase),
    stateful_step_count=50,
    report_multiple_bugs=True,
    suppress_health_check=(),
    deadline=duration(milliseconds=200),
    print_blob=False,
    backend="hypothesis",
)
settings.register_profile("default", default)
settings.load_profile("default")

assert settings.default is not None

CI = settings(
    derandomize=True,
    deadline=None,
    database=None,
    print_blob=True,
    suppress_health_check=[HealthCheck.too_slow],
)

settings.register_profile("ci", CI)


if is_in_ci():  # pragma: no cover # covered in ci, but not locally
    settings.load_profile("ci")

assert settings.default is not None


# Check that the kwonly args to settings.__init__ is the same as the set of
# defined settings - in case we've added or remove something from one but
# not the other.
assert set(all_settings) == {
    p.name
    for p in inspect.signature(settings.__init__).parameters.values()
    if p.kind == inspect.Parameter.KEYWORD_ONLY
}
