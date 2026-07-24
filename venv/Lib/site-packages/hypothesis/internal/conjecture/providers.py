# This file is part of Hypothesis, which may be found at
# https://github.com/HypothesisWorks/hypothesis/
#
# Copyright the Hypothesis Authors.
# Individual contributors are listed in AUTHORS.rst and the git log.
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at https://mozilla.org/MPL/2.0/.

import abc
import contextlib
import math
import sys
import warnings
from collections.abc import Iterable
from contextlib import AbstractContextManager, contextmanager
from functools import cached_property
from random import Random
from sys import float_info
from typing import (
    TYPE_CHECKING,
    Any,
    ClassVar,
    Literal,
    Optional,
    TypeAlias,
    TypedDict,
    TypeVar,
)

from sortedcontainers import SortedSet

from hypothesis.errors import HypothesisWarning
from hypothesis.internal.cache import LRUCache
from hypothesis.internal.compat import WINDOWS, int_from_bytes
from hypothesis.internal.conjecture.choice import (
    ChoiceConstraintsT,
    ChoiceT,
    ChoiceTypeT,
    FloatConstraints,
    choice_constraints_key,
    choice_permitted,
)
from hypothesis.internal.conjecture.floats import lex_to_float
from hypothesis.internal.conjecture.junkdrawer import bits_to_bytes
from hypothesis.internal.conjecture.utils import (
    Sampler,
    many,
)
from hypothesis.internal.constants_ast import (
    Constants,
    constants_from_module,
    is_local_module_file,
)
from hypothesis.internal.floats import (
    SIGNALING_NAN,
    clamp,
    float_to_int,
    make_float_clamper,
    next_down,
    next_up,
)
from hypothesis.internal.intervalsets import IntervalSet
from hypothesis.internal.observability import InfoObservationType, TestCaseObservation
from hypothesis.internal.statistics import (
    LogStudentTDistribution,
    PiecewiseDistribution,
    UniformDistribution,
)

if TYPE_CHECKING:
    from hypothesis.internal.conjecture.data import ConjectureData
    from hypothesis.internal.constants_ast import ConstantT

T = TypeVar("T")
LifetimeT: TypeAlias = Literal["test_case", "test_function"]
COLLECTION_DEFAULT_MAX_SIZE = 10**10  # "arbitrarily large"

# see https://github.com/HypothesisWorks/hypothesis/pull/4728 for details
INTEGERS_DISTRIBUTION = PiecewiseDistribution(
    inner=UniformDistribution(half_width=256),
    outer=LogStudentTDistribution(scale_bits=13, df=2),
    switchover=256,
)


#: Registered Hypothesis backends. This is a dictionary where keys are the name
#: to be used in |settings.backend|. The value of a key can be either:
#:
#: * A string corresponding to an importable absolute path of a
#:   |PrimitiveProvider| subclass
#: * A |PrimitiveProvider| subclass (the class itself, not an instance of the
#:   class)
#:
#: Hypothesis will instantiate the corresponding |PrimitiveProvider| subclass
#: when the backend is requested by a test's |settings.backend| value.
#:
#: For example, the default Hypothesis backend is registered as:
#:
#: .. code-block:: python
#:
#:    from hypothesis.internal.conjecture.providers import AVAILABLE_PROVIDERS
#:
#:    AVAILABLE_PROVIDERS["hypothesis"] = "hypothesis.internal.conjecture.providers.HypothesisProvider"
#:    # or
#:    AVAILABLE_PROVIDERS["hypothesis"] = HypothesisProvider
#:
#: And can be used with:
#:
#: .. code-block:: python
#:
#:     from hypothesis import given, settings, strategies as st
#:
#:     @given(st.integers())
#:     @settings(backend="hypothesis")
#:     def f(n):
#:         pass
#:
#: Though, as ``backend="hypothesis"`` is the default setting, the above would
#: typically not have any effect.
#:
#: For third-party backend authors, we strongly encourage ensuring that
#: ``import hypothesis`` does not automatically import the expensive parts of
#: your package, by:
#:
#: - setting a string path here, instead of a provider class
#: - ensuring the registered hypothesis plugin path references a path which just
#:   sets AVAILABLE_PROVIDERS and does not import your package
AVAILABLE_PROVIDERS: dict[str, str | type["PrimitiveProvider"]] = {
    "hypothesis": "hypothesis.internal.conjecture.providers.HypothesisProvider",
    "hypothesis-urandom": "hypothesis.internal.conjecture.providers.URandomProvider",
}
# cache the choice_permitted constants for a particular set of constraints.
CacheKeyT: TypeAlias = tuple[ChoiceTypeT, tuple[Any, ...]]
CacheValueT: TypeAlias = tuple[tuple["ConstantT", ...], tuple["ConstantT", ...]]
CONSTANTS_CACHE: LRUCache[CacheKeyT, CacheValueT] = LRUCache(1024)


_constants_integers = (
    # powers of 2
    [2**n for n in range(16, 66)]
    # powers of 10
    + [10**n for n in range(5, 20)]
    # factorials
    + [math.factorial(n) for n in range(9, 21)]
    # a few primorial numbers https://en.wikipedia.org/wiki/Primorial
    + [
        510510,
        6469693230,
        304250263527210,
        32589158477190044730,
    ]
)
_constants_integers.extend(
    [n - 1 for n in _constants_integers] + [n + 1 for n in _constants_integers]
)
_constants_integers.extend([-x for x in _constants_integers])

# arbitrary cutoffs to keep our list bounded
assert all(50_000 <= abs(n) <= 2**66 for n in _constants_integers)

_constant_floats = (
    [
        0.5,
        1.1,
        1.5,
        1.9,
        1.0 / 3,
        10e6,
        10e-6,
        1.175494351e-38,
        next_up(0.0),
        float_info.min,
        float_info.max,
        3.402823466e38,
        9007199254740992.0,
        1 - 10e-6,
        2 + 10e-6,
        1.192092896e-07,
        2.2204460492503131e-016,
    ]
    + [2.0**-n for n in (24, 14, 149, 126)]  # minimum (sub)normals for float16,32
    + [float_info.min / n for n in (2, 10, 1000, 100_000)]  # subnormal in float64
)
_constant_floats.extend([-x for x in _constant_floats])
assert all(isinstance(f, float) for f in _constant_floats)

_constant_strings = {
    # strings which can be interpreted as code / logic
    "undefined",
    "null",
    "NULL",
    "nil",
    "NIL",
    "true",
    "false",
    "True",
    "False",
    "TRUE",
    "FALSE",
    "None",
    "none",
    "if",
    "then",
    "else",
    "__dict__",
    "__proto__",  # javascript
    # strings which can be interpreted as a number
    "0",
    "1e100",
    "0..0",
    "0/0",
    "1/0",
    "+0.0",
    "Infinity",
    "-Infinity",
    "Inf",
    "INF",
    "NaN",
    "9" * 30,
    # common ascii characters
    ",./;'[]\\-=<>?:\"{}|_+!@#$%^&*()`~",
    # common unicode characters
    "Ω≈ç√∫˜µ≤≥÷åß∂ƒ©˙∆˚¬…æœ∑´®†¥¨ˆøπ“‘¡™£¢∞§¶•ªº–≠¸˛Ç◊ı˜Â¯˘¿ÅÍÎÏ˝ÓÔÒÚÆ☃Œ„´‰ˇÁ¨ˆØ∏”’`⁄€‹›ﬁﬂ‡°·‚—±",
    # characters which increase in length when lowercased
    "Ⱥ",
    "Ⱦ",
    # ligatures
    "æœÆŒﬀʤʨß",
    # emoticons
    "(╯°□°）╯︵ ┻━┻)",
    # emojis
    "😍",
    "🇺🇸",
    # emoji modifiers
    "🏻",  # U+1F3FB Light Skin Tone,
    "👍🏻",  # 👍 followed by U+1F3FB
    # RTL text
    "الكل في المجمو عة",
    # Ogham text, which contains the only character in the Space Separators
    # unicode category (Zs) that isn't visually blank:  .  # noqa: RUF003
    "᚛ᚄᚓᚐᚋᚒᚄ ᚑᚄᚂᚑᚏᚅ᚜",
    # thai consonant + spacing vowel combinations, which have unusual visual combining behavior
    "กา",
    "ก ำกำ",
    # readable variations on text (bolt/italic/script)
    "𝐓𝐡𝐞 𝐪𝐮𝐢𝐜𝐤 𝐛𝐫𝐨𝐰𝐧 𝐟𝐨𝐱 𝐣𝐮𝐦𝐩𝐬 𝐨𝐯𝐞𝐫 𝐭𝐡𝐞 𝐥𝐚𝐳𝐲 𝐝𝐨𝐠",
    "𝕿𝖍𝖊 𝖖𝖚𝖎𝖈𝖐 𝖇𝖗𝖔𝖜𝖓 𝖋𝖔𝖝 𝖏𝖚𝖒𝖕𝖘 𝖔𝖛𝖊𝖗 𝖙𝖍𝖊 𝖑𝖆𝖟𝖞 𝖉𝖔𝖌",
    "𝑻𝒉𝒆 𝒒𝒖𝒊𝒄𝒌 𝒃𝒓𝒐𝒘𝒏 𝒇𝒐𝒙 𝒋𝒖𝒎𝒑𝒔 𝒐𝒗𝒆𝒓 𝒕𝒉𝒆 𝒍𝒂𝒛𝒚 𝒅𝒐𝒈",
    "𝓣𝓱𝓮 𝓺𝓾𝓲𝓬𝓴 𝓫𝓻𝓸𝔀𝓷 𝓯𝓸𝔁 𝓳𝓾𝓶𝓹𝓼 𝓸𝓿𝓮𝓻 𝓽𝓱𝓮 𝓵𝓪𝔃𝔂 𝓭𝓸𝓰",
    "𝕋𝕙𝕖 𝕢𝕦𝕚𝕔𝕜 𝕓𝕣𝕠𝕨𝕟 𝕗𝕠𝕩 𝕛𝕦𝕞𝕡𝕤 𝕠𝕧𝕖𝕣 𝕥𝕙𝕖 𝕝𝕒𝕫𝕪 𝕕𝕠𝕘",
    # upsidown text
    "ʇǝɯɐ ʇᴉs ɹolop ɯnsdᴉ ɯǝɹo˥",
    # reserved strings in windows
    "NUL",
    "COM1",
    "LPT1",
    # scunthorpe problem
    "Scunthorpe",
    # zalgo text
    "Ṱ̺̺̕o͞ ̷i̲̬͇̪͙n̝̗͕v̟̜̘̦͟o̶̙̰̠kè͚̮̺̪̹̱̤ ̖t̝͕̳̣̻̪͞h̼͓̲̦̳̘̲e͇̣̰̦̬͎ ̢̼̻̱̘h͚͎͙̜̣̲ͅi̦̲̣̰̤v̻͍e̺̭̳̪̰-m̢iͅn̖̺̞̲̯̰d̵̼̟͙̩̼̘̳ ̞̥̱̳̭r̛̗̘e͙p͠r̼̞̻̭̗e̺̠̣͟s̘͇̳͍̝͉e͉̥̯̞̲͚̬͜ǹ̬͎͎̟̖͇̤t͍̬̤͓̼̭͘ͅi̪̱n͠g̴͉ ͏͉ͅc̬̟h͡a̫̻̯͘o̫̟̖͍̙̝͉s̗̦̲.̨̹͈̣",
    #
    # examples from https://faultlore.com/blah/text-hates-you/
    "मनीष منش",
    "पन्ह पन्ह त्र र्च कृकृ ड्ड न्हृे إلا بسم الله",
    "lorem لا بسم الله ipsum 你好1234你好",
    # unicode charaacters causing unconditional line breaks, as defined by UAX #14:
    # https://www.unicode.org/reports/tr14/.
    #
    # We've seen multiple bugs caused by assuming `str.splitlines` is equivalent to
    # splitting over "\n", while it actually splits over all line breaks!
    #
    # We intersperse the line breaks with normal characters to increase the likelihood
    # of triggering such a bug.
    (
        "a"
        "\u000a"  # line feed (class: LF)
        "b"
        "\u000d"  # carriage return (class: CR)
        "c"
        "\u0085"  # next line (class: NL)
        "d"
        "\u000b"  # line tabulation (class: BK)
        "e"
        "\u000c"  # form feed (class: BK)
        "f"
        "\u2028"  # line separator (class: BK)
        "g"
        "\u2029"  # paragraph separator (class: BK)
        "h"
        "\u000d\u000a"  # CR+LF
        "i"
    ),
}


# we don't actually care what order the constants are sorted in, just that the
# ordering is deterministic.
GLOBAL_CONSTANTS = Constants(
    integers=SortedSet(_constants_integers),
    floats=SortedSet(_constant_floats, key=float_to_int),
    bytes=SortedSet(),
    strings=SortedSet(_constant_strings),
)

_local_constants = Constants(
    integers=SortedSet(),
    floats=SortedSet(key=float_to_int),
    bytes=SortedSet(),
    strings=SortedSet(),
)
# Modules that we've already seen and processed for local constants. These are
# are all modules, not necessarily local ones. This lets us quickly see which
# modules are new without an expensive path.resolve() or is_local_module_file
# cache lookup.
#
# We track by id so we can handle users that put unhashable types like SimpleNamespace
# into sys.modules. ModuleType.__hash__ falls back to id, so this is equivalent
# in the standard case.
_seen_modules: set[int] = set()
_sys_modules_len: int | None = None


def _get_local_constants() -> Constants:
    global _sys_modules_len, _local_constants

    if sys.platform == "emscripten":  # pragma: no cover
        # pyodide builds bundle the stdlib in a nonstandard location, like
        # `/lib/python312.zip/heapq.py`. To avoid identifying the entirety of
        # the stdlib as local code and slowing down on emscripten, instead return
        # that nothing is local.
        #
        # pyodide may provide some way to distinguish stdlib/third-party/local
        # code. I haven't looked into it. If they do, we should correctly implement
        # ModuleLocation for pyodide instead of this.
        return _local_constants

    count_constants = len(_local_constants)
    # We call this function once per HypothesisProvider instance, i.e. once per
    # input, so it needs to be performant. The logic here is more complicated
    # than necessary because of this.
    #
    # First, we check whether there are any new modules with a very cheap length
    # check. This check can be fooled if a module is added while another module is
    # removed, but the more correct check against tuple(sys.modules.keys()) is
    # substantially more expensive. Such a new module would eventually be discovered
    # if / when the length changes again in the future.
    #
    # If the length has changed, we find just modules we haven't seen before. Of
    # those, we find the ones which correspond to local modules, and extract their
    # constants.

    # careful: store sys.modules length when we first check to avoid race conditions
    # with other threads loading a module before we set _sys_modules_len.
    if (sys_modules_len := len(sys.modules)) != _sys_modules_len:
        new_modules = [
            module
            for module in list(sys.modules.values())
            if id(module) not in _seen_modules
        ]
        # Repeated SortedSet unions are expensive. Do the initial unions on a
        # set(), then do a one-time union with _local_constants after.
        new_constants = Constants()
        for module in new_modules:
            if (
                module_file := getattr(module, "__file__", None)
            ) is not None and is_local_module_file(module_file):
                new_constants |= constants_from_module(module)
            _seen_modules.add(id(module))
        _local_constants |= new_constants
        _sys_modules_len = sys_modules_len

    # if we add any new constant, invalidate the constant cache for permitted values.
    # A more efficient approach would be invalidating just the keys with this
    # choice_type.
    if len(_local_constants) > count_constants:
        CONSTANTS_CACHE.cache.clear()

    return _local_constants


@contextmanager
def with_register_backend(name, provider_cls):
    try:
        AVAILABLE_PROVIDERS[name] = provider_cls
        yield
    finally:
        del AVAILABLE_PROVIDERS[name]


class _BackendInfoMsg(TypedDict):
    type: InfoObservationType
    title: str
    content: str | dict[str, Any]


# TODO_DOCS: link to choice sequence explanation page


class PrimitiveProvider(abc.ABC):
    """
    |PrimitiveProvider| is the implementation interface of a
    :ref:`Hypothesis backend <alternative-backends>`.

    A |PrimitiveProvider| is required to implement the following five
    ``draw_*`` methods:

    * |PrimitiveProvider.draw_integer|
    * |PrimitiveProvider.draw_boolean|
    * |PrimitiveProvider.draw_float|
    * |PrimitiveProvider.draw_string|
    * |PrimitiveProvider.draw_bytes|

    Each strategy in Hypothesis generates values by drawing a series of choices
    from these five methods. By overriding them, a |PrimitiveProvider| can control
    the distribution of inputs generated by Hypothesis.

    For example, :pypi:`hypothesis-crosshair` implements a |PrimitiveProvider|
    which uses an SMT solver to generate inputs that uncover new branches.

    Once you implement a |PrimitiveProvider|, you can make it available for use
    through |AVAILABLE_PROVIDERS|.
    """

    #: The lifetime of a |PrimitiveProvider| instance. Either ``test_function``
    #: or ``test_case``.
    #:
    #: If ``test_function`` (the default), a single provider instance will be
    #: instantiated and used for the entirety of each test function (i.e., roughly
    #: one provider per |@given| annotation). This can be useful for tracking state
    #: over the entirety of a test function.
    #:
    #: If ``test_case``, a new provider instance will be instantiated and used for
    #: each input Hypothesis generates.
    #:
    #: The ``conjecturedata`` argument to ``PrimitiveProvider.__init__`` will
    #: be ``None`` for a lifetime of ``test_function``, and an instance of
    #: ``ConjectureData`` for a lifetime of ``test_case``.
    #:
    #: Third-party providers likely want to set a lifetime of ``test_function``.
    lifetime: ClassVar[LifetimeT] = "test_function"

    #: Solver-based backends such as ``hypothesis-crosshair`` use symbolic values
    #: which record operations performed on them in order to discover new paths.
    #: If ``avoid_realization`` is set to ``True``, hypothesis will avoid interacting
    #: with symbolic choices returned by the provider in any way that would force
    #: the solver to narrow the range of possible values for that symbolic.
    #:
    #: Setting this to ``True`` disables some hypothesis features and optimizations.
    #: Only set this to ``True`` if it is necessary for your backend.
    avoid_realization: ClassVar[bool] = False

    #: If ``True``, |PrimitiveProvider.on_observation| will be added as a
    #: callback via |add_observability_callback|, enabling observability during
    # the lifetime of this provider. If ``False``, |PrimitiveProvider.on_observation|
    #: will never be called by Hypothesis.
    #:
    #: The opt-in behavior of observability is because enabling observability
    #: might increase runtime or memory usage.
    add_observability_callback: ClassVar[bool] = False

    def __init__(self, conjecturedata: Optional["ConjectureData"], /) -> None:
        self._cd = conjecturedata

    @abc.abstractmethod
    def draw_boolean(
        self,
        p: float = 0.5,
    ) -> bool:
        """
        Draw a boolean choice.

        Parameters
        ----------
        p: float
            The probability of returning ``True``. Between 0 and 1 inclusive.

            Except for ``0`` and ``1``, the value of ``p`` is a hint provided by
            Hypothesis, and may be ignored by the backend.

            If ``0``, the provider must return ``False``. If ``1``, the provider
            must return ``True``.
        """
        raise NotImplementedError

    @abc.abstractmethod
    def draw_integer(
        self,
        min_value: int | None = None,
        max_value: int | None = None,
        *,
        weights: dict[int, float] | None = None,
        shrink_towards: int = 0,
    ) -> int:
        """
        Draw an integer choice.

        Parameters
        ----------
        min_value : int | None
            (Inclusive) lower bound on the integer value. If ``None``, there is
            no lower bound.
        max_value : int | None
            (Inclusive) upper bound on the integer value. If ``None``, there is
            no upper bound.
        weights: dict[int, float] | None
            Maps keys in the range [``min_value``, ``max_value``] to the probability
            of returning that key.
        shrink_towards: int
            The integer to shrink towards. This is not used during generation and
            can be ignored by backends.
        """
        raise NotImplementedError

    @abc.abstractmethod
    def draw_float(
        self,
        *,
        min_value: float = -math.inf,
        max_value: float = math.inf,
        allow_nan: bool = True,
        smallest_nonzero_magnitude: float,
    ) -> float:
        """
        Draw a float choice.

        Parameters
        ----------
        min_value : float
            (Inclusive) lower bound on the float value.
        max_value : float
            (Inclusive) upper bound on the float value.
        allow_nan : bool
            If ``False``, it is invalid to return ``math.nan``.
        smallest_nonzero_magnitude : float
            The smallest allowed nonzero magnitude. ``draw_float`` should not
            return a float ``f`` if ``abs(f) < smallest_nonzero_magnitude``.
        """
        raise NotImplementedError

    @abc.abstractmethod
    def draw_string(
        self,
        intervals: IntervalSet,
        *,
        min_size: int = 0,
        max_size: int = COLLECTION_DEFAULT_MAX_SIZE,
    ) -> str:
        """
        Draw a string choice.

        Parameters
        ----------
        intervals : IntervalSet
            The set of codepoints to sample from.
        min_size : int
            (Inclusive) lower bound on the string length.
        max_size : int
            (Inclusive) upper bound on the string length.
        """
        raise NotImplementedError

    @abc.abstractmethod
    def draw_bytes(
        self,
        min_size: int = 0,
        max_size: int = COLLECTION_DEFAULT_MAX_SIZE,
    ) -> bytes:
        """
        Draw a bytes choice.

        Parameters
        ----------
        min_size : int
            (Inclusive) lower bound on the bytes length.
        max_size : int
            (Inclusive) upper bound on the bytes length.
        """
        raise NotImplementedError

    def per_test_case_context_manager(self) -> AbstractContextManager:
        """
        Returns a context manager which will be entered each time Hypothesis
        starts generating and executing one |test case|, and exited when that test
        case finishes generating and executing, including if any exception is
        thrown.

        In the lifecycle of a Hypothesis test, this is called before
        generating strategy values for each test case. This is just before any
        :ref:`custom executor <custom-function-execution>` is called.

        Even if not returning a custom context manager, |PrimitiveProvider|
        subclasses are welcome to override this method to know when Hypothesis
        starts and ends the execution of a single test case.
        """
        return contextlib.nullcontext()

    def realize(self, value: T, *, for_failure: bool = False) -> T:
        """
        Called whenever hypothesis requires a concrete (non-symbolic) value from
        a potentially symbolic value. Hypothesis will not check that ``value`` is
        symbolic before calling ``realize``, so you should handle the case where
        ``value`` is non-symbolic.

        The returned value should be non-symbolic.  If you cannot provide a value,
        raise |BackendCannotProceed| with a value of ``"discard_test_case"``.

        If ``for_failure`` is ``True``, the value is associated with a |failing test case|.
        In this case, the backend should spend substantially more effort when
        attempting to realize the value, since it is important to avoid discarding
        failing test cases. Backends may still raise |BackendCannotProceed| when
        ``for_failure`` is ``True``, if realization is truly impossible or if
        realization takes significantly longer than expected (say, 5 minutes).
        """
        return value

    def replay_choices(self, choices: tuple[ChoiceT, ...]) -> None:
        """
        Called when Hypothesis has discovered a |choice sequence| which the provider
        may wish to enqueue to replay under its own instrumentation when we next
        ask to generate a |test case|, rather than generating one from scratch.

        This is used to e.g. warm-start :pypi:`hypothesis-crosshair` with a corpus
        of high-code-coverage inputs discovered by
        `HypoFuzz <https://hypofuzz.com/>`_.
        """
        return None

    def observe_test_case(self) -> dict[str, Any]:
        """Called at the end of the |test case| when :ref:`observability
        <observability>` is enabled.

        The return value should be a non-symbolic json-encodable dictionary,
        and will be included in observations as ``observation["metadata"]["backend"]``.
        """
        return {}

    def observe_information_messages(
        self, *, lifetime: LifetimeT
    ) -> Iterable[_BackendInfoMsg]:
        """Called at the end of each |test case| and again at end of the test function.

        Return an iterable of ``{type: info/alert/error, title: str, content: str | dict}``
        dictionaries to be delivered as individual information messages. Hypothesis
        adds the ``run_start`` timestamp and ``property`` name for you.
        """
        assert lifetime in ("test_case", "test_function")
        yield from []

    def on_observation(self, observation: TestCaseObservation) -> None:  # noqa: B027
        """
        Called at the end of each |test case| which uses this provider, with the same
        ``observation["type"] == "test_case"`` observation that is passed to
        other callbacks added via |add_observability_callback|. This method is not
        called with ``observation["type"] in {"info", "alert", "error"}``
        observations.

        .. important::

            For |PrimitiveProvider.on_observation| to be called by Hypothesis,
            |PrimitiveProvider.add_observability_callback| must be set to ``True``.

            |PrimitiveProvider.on_observation| is explicitly opt-in, as enabling
            observability might increase runtime or memory usage.

        Calls to this method are guaranteed to alternate with calls to
        |PrimitiveProvider.per_test_case_context_manager|. For example:

        .. code-block:: python

            # test function starts
            per_test_case_context_manager()
            on_observation()
            per_test_case_context_manager()
            on_observation()
            ...
            # test function ends

        Note that |PrimitiveProvider.on_observation| will not be called for test
        cases which did not use this provider during generation, for example
        during |Phase.reuse| or |Phase.shrink|, or because Hypothesis switched
        to the standard Hypothesis backend after this backend raised too many
        |BackendCannotProceed| exceptions.
        """

    def span_start(self, label: int, /) -> None:  # noqa: B027  # non-abstract noop
        """Marks the beginning of a semantically meaningful span of choices.

        Spans are a depth-first tree structure. A span is opened by a call to
        |PrimitiveProvider.span_start|, and a call to |PrimitiveProvider.span_end|
        closes the most recently opened span. So the following sequence of calls:

        .. code-block:: python

            span_start(label=1)
            n1 = draw_integer()
            span_start(label=2)
            b1 = draw_boolean()
            n2 = draw_integer()
            span_end()
            f1 = draw_float()
            span_end()

        produces the following two spans of choices:

        .. code-block::

            1: [n1, b1, n2, f1]
            2: [b1, n2]

        Hypothesis uses spans to denote "semantically meaningful" sequences of
        choices. For instance, Hypothesis opens a span for the sequence of choices
        made while drawing from each strategy. Not every span corresponds to a
        strategy; the generation of e.g. each element in |st.lists| is also marked
        with a span, among others.

        ``label`` is an opaque integer, which has no defined semantics.
        The only guarantee made by Hypothesis is that all spans with the same
        "meaning" will share the same ``label``. So all spans from the same
        strategy will share the same label, as will e.g. the spans for |st.lists|
        elements.

        Providers can track calls to |PrimitiveProvider.span_start| and
        |PrimitiveProvider.span_end| to learn something about the semantics of
        the test's |choice sequence|. For instance, a provider could track the depth
        of the span tree, or the number of unique labels, which says something about
        the complexity of the choices being generated. Or a provider could track
        the span tree across |test cases| in order to determine what strategies are
        being used in what contexts.

        It is possible for Hypothesis to start and immediately stop a span,
        without calling a ``draw_*`` method in between. These spans contain zero
        choices.

        Hypothesis will always balance the number of calls to
        |PrimitiveProvider.span_start| and |PrimitiveProvider.span_end|. A call
        to |PrimitiveProvider.span_start| will always be followed by a call to
        |PrimitiveProvider.span_end| before the end of the test case.

        |PrimitiveProvider.span_start| is called from ``ConjectureData.start_span()``
        internally.
        """

    def span_end(self, discard: bool, /) -> None:  # noqa: B027
        """Marks the end of a semantically meaningful span of choices.

        ``discard`` is ``True`` when the draw was filtered out or otherwise marked
        as unlikely to contribute to the input data as seen by the user's test.
        Note however that side effects can make this determination unsound.

        |PrimitiveProvider.span_end| is called from ``ConjectureData.stop_span()``
        internally.
        """


class HypothesisProvider(PrimitiveProvider):
    lifetime = "test_case"

    def __init__(self, conjecturedata: Optional["ConjectureData"], /):
        super().__init__(conjecturedata)
        self._random = None if self._cd is None else self._cd._random

    @cached_property
    def _local_constants(self):
        # defer computation of local constants until/if we need it
        return _get_local_constants()

    def _maybe_draw_constant(
        self,
        choice_type: ChoiceTypeT,
        constraints: ChoiceConstraintsT,
        *,
        p: float = 0.05,
    ) -> Optional["ConstantT"]:
        assert self._random is not None
        assert choice_type != "boolean"
        # check whether we even want a constant before spending time computing
        # and caching the allowed constants.
        if self._random.random() > p:
            return None

        # note: this property access results in computation being done
        assert self._local_constants is not None

        key = (choice_type, choice_constraints_key(choice_type, constraints))
        if key not in CONSTANTS_CACHE:
            CONSTANTS_CACHE[key] = (
                tuple(
                    choice
                    for choice in GLOBAL_CONSTANTS.set_for_type(choice_type)
                    if choice_permitted(choice, constraints)
                ),
                tuple(
                    choice
                    for choice in self._local_constants.set_for_type(choice_type)
                    if choice_permitted(choice, constraints)
                ),
            )

        # split constants into two pools, so we still have a good chance to draw
        # global constants even if there are many local constants.
        global_constants, local_constants = CONSTANTS_CACHE[key]
        constants_lists = ([global_constants] if global_constants else []) + (
            [local_constants] if local_constants else []
        )
        if not constants_lists:
            return None

        # At this point, we've decided to use a constant. Now we select which pool
        # to draw that constant from.
        #
        # Note that this approach has a different probability distribution than
        # attempting a random.random for both global_constants and local_constants.
        constants = self._random.choice(constants_lists)
        return self._random.choice(constants)

    def draw_boolean(
        self,
        p: float = 0.5,
    ) -> bool:
        assert self._random is not None

        if p <= 0:
            return False
        if p >= 1:
            return True

        return self._random.random() < p

    def draw_integer(
        self,
        min_value: int | None = None,
        max_value: int | None = None,
        *,
        weights: dict[int, float] | None = None,
        shrink_towards: int = 0,
    ) -> int:
        assert self._cd is not None
        if (
            constant := self._maybe_draw_constant(
                "integer",
                {
                    "min_value": min_value,
                    "max_value": max_value,
                    "weights": weights,
                    "shrink_towards": shrink_towards,
                },
            )
        ) is not None:
            assert isinstance(constant, int)
            return constant

        if weights is not None:
            assert min_value is not None
            assert max_value is not None

            # format of weights is a mapping of ints to p, where sum(p) < 1.
            # The remaining probability mass is uniformly distributed over
            # *all* ints (not just the unmapped ones; this is somewhat undesirable,
            # but simplifies things).
            #
            # We assert that sum(p) is strictly less than 1 because it simplifies
            # handling forced values when we can force into the unmapped probability
            # mass. We should eventually remove this restriction.
            sampler = Sampler(
                [1 - sum(weights.values()), *weights.values()], observe=False
            )
            # if we're forcing, it's easiest to force into the unmapped probability
            # mass and then force the drawn value after.
            idx = sampler.sample(self._cd)

            if idx == 0:
                return self._draw_integer_from_distribution(min_value, max_value)
            # implicit reliance on dicts being sorted for determinism
            return list(weights)[idx - 1]

        return self._draw_integer_from_distribution(min_value, max_value)

    def _draw_integer_from_distribution(
        self, min_value: int | None, max_value: int | None
    ) -> int:
        assert self._random is not None
        dist = INTEGERS_DISTRIBUTION

        # Our integers distribution is defined over the full float64 range. If the user
        # has not placed a lower and/or upper bound, we don't actually want to sample
        # from this full range, because we might otherwise generate integers
        # that are so large as to cause problems. For example python errors when passing
        # a too-large integer to c APIs (common and implicit in some python stdlib APIs).
        #
        # * If the user asks for an unbounded range, we bound at [-2**128, 2**128].
        # * If the user asks for a half-bounded range starting at n, we bound at
        #   [n, 2**128] for n < 2**127, and [n, 2n] otherwise, so we always support
        #   generating large integers if a user explicitly asks for integers around that
        #   magnitude. The choice of a factor of two here is arbitrary.
        if min_value is None and max_value is None:
            min_value = -(2**128)
            max_value = 2**128
        elif min_value is None:
            assert max_value is not None
            min_value = -max(2**128, 2 * abs(max_value))
        elif max_value is None:
            max_value = max(2**128, 2 * abs(min_value))

        lo = 0.0
        hi = 1.0
        safe_bounds = True
        try:
            cdf_min = min_value - 0.5
            cdf_max = max_value + 0.5
        except OverflowError:
            safe_bounds = False

        if safe_bounds:
            lo = dist.cdf(cdf_min)
            hi = dist.cdf(cdf_max)

            # if the bounds in CDF-space are too close together, there isn't enough room
            # to stably sample between hi and lo. For example, in the extreme case of
            # hi = lo (which can happen even if min_value != max_value due to either float
            # imprecision or numerical instability), our sampling algorithm would collapse
            # to always return the sample value.
            #
            # The choice of 1e-13 here is somewhat arbitrary. At 1.0, epsilon in float64
            # is 2^-53 ~= 1e-16. We want some breathing room from epsilon, because we need
            # some variation to actually sample in. I haven't put much thought into where
            # the correct tradeoff lies. The downside risk to choosing a too-aggressive
            # bound is some integers may literally not be representable in the `hi - lo`
            # sampling space. The downside risk to a too-conservative bound is switching
            # off our good distribution too early.
            if hi - lo < 1e-13:
                safe_bounds = False

        if safe_bounds:
            # inverse_cdf requires strictly 0 < p < 1. Resample until we get it.
            while (p := lo + self._random.random() * (hi - lo)) in {0, 1}:
                pass  # pragma: no cover
            n = round(dist.inverse_cdf(p))
            # As an extra measure of safety, clamp our generated value to the requested
            # range. It would of course be nice if we provably satisfied this by
            # construction - I'm just not 100% confident that we do.
            return clamp(min_value, n, max_value)

        # We have determined the bounds [min_value, max_value] are not numerically safe
        # to use. Fall back to a simpler and safer distribution: a uniform distribution
        # on [min_value, max_value].
        return self._random.randint(min_value, max_value)

    def draw_float(
        self,
        *,
        min_value: float = -math.inf,
        max_value: float = math.inf,
        allow_nan: bool = True,
        smallest_nonzero_magnitude: float,
    ) -> float:
        assert self._random is not None

        constraints: FloatConstraints = {
            "min_value": min_value,
            "max_value": max_value,
            "allow_nan": allow_nan,
            "smallest_nonzero_magnitude": smallest_nonzero_magnitude,
        }
        if (
            constant := self._maybe_draw_constant("float", constraints, p=0.15)
        ) is not None:
            assert isinstance(constant, float)
            return constant

        # on top of the probability to draw a constant float, we independently
        # upweight 0.0/-0.0, math.inf, -math.inf, nans, and boundary values.
        weird_floats = [
            f
            for f in [
                0.0,
                -0.0,
                math.inf,
                -math.inf,
                math.nan,
                -math.nan,
                SIGNALING_NAN,
                -SIGNALING_NAN,
                min_value,
                next_up(min_value),
                min_value + 1,
                max_value - 1,
                next_down(max_value),
                max_value,
            ]
            if choice_permitted(f, constraints)
        ]

        if weird_floats and self._random.random() < 0.05:
            return self._random.choice(weird_floats)

        clamper = make_float_clamper(
            min_value,
            max_value,
            smallest_nonzero_magnitude=smallest_nonzero_magnitude,
            allow_nan=allow_nan,
        )

        result = self._draw_float()
        if allow_nan and math.isnan(result):
            clamped = result  # pragma: no cover
        else:
            clamped = clamper(result)
        if float_to_int(clamped) != float_to_int(result) and not (
            math.isnan(result) and allow_nan
        ):
            result = clamped
        return result

    def draw_string(
        self,
        intervals: IntervalSet,
        *,
        min_size: int = 0,
        max_size: int = COLLECTION_DEFAULT_MAX_SIZE,
    ) -> str:
        assert self._cd is not None
        assert self._random is not None

        if len(intervals) == 0:
            return ""

        if (
            constant := self._maybe_draw_constant(
                "string",
                {"intervals": intervals, "min_size": min_size, "max_size": max_size},
            )
        ) is not None:
            assert isinstance(constant, str)
            return constant

        average_size = min(
            max(min_size * 2, min_size + 5),
            0.5 * (min_size + max_size),
        )

        chars = []
        elements = many(
            self._cd,
            min_size=min_size,
            max_size=max_size,
            average_size=average_size,
            observe=False,
        )
        while elements.more():
            if len(intervals) > 256:
                if self.draw_boolean(0.2):
                    i = self._random.randint(256, len(intervals) - 1)
                else:
                    i = self._random.randint(0, 255)
            else:
                i = self._random.randint(0, len(intervals) - 1)

            chars.append(intervals.char_in_shrink_order(i))

        return "".join(chars)

    def draw_bytes(
        self,
        min_size: int = 0,
        max_size: int = COLLECTION_DEFAULT_MAX_SIZE,
    ) -> bytes:
        assert self._cd is not None
        assert self._random is not None

        if (
            constant := self._maybe_draw_constant(
                "bytes", {"min_size": min_size, "max_size": max_size}
            )
        ) is not None:
            assert isinstance(constant, bytes)
            return constant

        buf = bytearray()
        average_size = min(
            max(min_size * 2, min_size + 5),
            0.5 * (min_size + max_size),
        )
        elements = many(
            self._cd,
            min_size=min_size,
            max_size=max_size,
            average_size=average_size,
            observe=False,
        )
        while elements.more():
            buf += self._random.randbytes(1)

        return bytes(buf)

    def _draw_float(self) -> float:
        assert self._random is not None

        f = lex_to_float(self._random.getrandbits(64))
        sign = 1 if self._random.getrandbits(1) else -1
        return sign * f


# Masks for masking off the first byte of an n-bit buffer.
# The appropriate mask is stored at position n % 8.
BYTE_MASKS = [(1 << n) - 1 for n in range(8)]
BYTE_MASKS[0] = 255


class BytestringProvider(PrimitiveProvider):
    lifetime = "test_case"

    def __init__(
        self, conjecturedata: Optional["ConjectureData"], /, *, bytestring: bytes
    ):
        super().__init__(conjecturedata)
        self.bytestring = bytestring
        self.index = 0
        self.drawn = bytearray()

    def _draw_bits(self, n):
        if n == 0:  # pragma: no cover
            return 0
        n_bytes = bits_to_bytes(n)
        if self.index + n_bytes > len(self.bytestring):
            self._cd.mark_overrun()
        buf = bytearray(self.bytestring[self.index : self.index + n_bytes])
        self.index += n_bytes

        buf[0] &= BYTE_MASKS[n % 8]
        buf = bytes(buf)
        self.drawn += buf
        return int_from_bytes(buf)

    def draw_boolean(
        self,
        p: float = 0.5,
    ) -> bool:
        if p <= 0:
            return False
        if p >= 1:
            return True

        # always use one byte for booleans to maintain constant draw size.
        # If a probability requires more than 8 bits to represent precisely,
        # the result will be slightly biased, but not badly.
        bits = 8
        size = 2**bits
        # always leave at least one value that can be true, even for very small
        # p.
        falsey = max(1, math.floor(size * (1 - p)))
        n = self._draw_bits(bits)
        return n >= falsey

    def draw_integer(
        self,
        min_value: int | None = None,
        max_value: int | None = None,
        *,
        weights: dict[int, float] | None = None,
        shrink_towards: int = 0,
    ) -> int:
        assert self._cd is not None

        # we explicitly ignore integer weights for now, as they are likely net
        # negative on fuzzer performance.

        if min_value is None and max_value is None:
            min_value = -(2**127)
            max_value = 2**127 - 1
        elif min_value is None:
            assert max_value is not None
            min_value = max_value - 2**64
        elif max_value is None:
            assert min_value is not None
            max_value = min_value + 2**64

        if min_value == max_value:
            return min_value

        bits = (max_value - min_value).bit_length()
        value = self._draw_bits(bits)
        while not (min_value <= value <= max_value):
            value = self._draw_bits(bits)
        return value

    def draw_float(
        self,
        *,
        min_value: float = -math.inf,
        max_value: float = math.inf,
        allow_nan: bool = True,
        smallest_nonzero_magnitude: float,
    ) -> float:
        n = self._draw_bits(64)
        sign = -1 if n >> 64 else 1
        f = sign * lex_to_float(n & ((1 << 64) - 1))
        clamper = make_float_clamper(
            min_value,
            max_value,
            smallest_nonzero_magnitude=smallest_nonzero_magnitude,
            allow_nan=allow_nan,
        )
        return clamper(f)

    def _draw_collection(self, min_size, max_size, *, alphabet_size):
        average_size = min(
            max(min_size * 2, min_size + 5),
            0.5 * (min_size + max_size),
        )
        elements = many(
            self._cd,
            min_size=min_size,
            max_size=max_size,
            average_size=average_size,
            observe=False,
        )
        values = []
        while elements.more():
            values.append(self.draw_integer(0, alphabet_size - 1))
        return values

    def draw_string(
        self,
        intervals: IntervalSet,
        *,
        min_size: int = 0,
        max_size: int = COLLECTION_DEFAULT_MAX_SIZE,
    ) -> str:
        values = self._draw_collection(min_size, max_size, alphabet_size=len(intervals))
        return "".join(chr(intervals[v]) for v in values)

    def draw_bytes(
        self,
        min_size: int = 0,
        max_size: int = COLLECTION_DEFAULT_MAX_SIZE,
    ) -> bytes:
        values = self._draw_collection(min_size, max_size, alphabet_size=2**8)
        return bytes(values)


class URandom(Random):
    # we reimplement a Random instance instead of using SystemRandom, because
    # os.urandom is not guaranteed to read from /dev/urandom.

    @staticmethod
    def _urandom(size: int) -> bytes:
        # By default, we buffer more data from /dev/urandom than the actual `size` number
        # of bytes we requested. This trips up anyone (and particularly Antithesis)
        # hooking /dev/urandom reads to control randomness.
        #
        # buffering=0 disables this behavior.
        with open("/dev/urandom", "rb", buffering=0) as f:
            return f.read(size)

    def getrandbits(self, k: int) -> int:
        assert k >= 0
        size = bits_to_bytes(k)
        n = int_from_bytes(self._urandom(size))
        # trim excess bits
        return n >> (size * 8 - k)

    def random(self) -> float:
        # adapted from random.SystemRandom.random
        return (int_from_bytes(self._urandom(7)) >> 3) * (2**-53)


class URandomProvider(HypothesisProvider):
    # A provider which reads directly from /dev/urandom as its source of randomness.
    # This provider exists to provide better Hypothesis integration with Antithesis
    # (https://antithesis.com/), which interprets calls to /dev/urandom as the
    # randomness to mutate. This effectively gives Antithesis control over
    # the choices made by the URandomProvider.
    #
    # If you are not using Antithesis, you probably don't want to use this
    # provider.

    def __init__(self, conjecturedata: Optional["ConjectureData"], /):
        super().__init__(conjecturedata)
        if WINDOWS:  # pragma: no cover
            warnings.warn(
                "/dev/urandom is not available on windows. Falling back to "
                'standard PRNG generation (equivalent to backend="hypothesis").',
                HypothesisWarning,
                stacklevel=1,
            )
            # don't overwrite the HypothesisProvider self._random attribute in
            # this case
        else:
            self._random = URandom()
