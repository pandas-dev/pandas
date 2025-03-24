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
import warnings
from collections.abc import Iterable
from random import Random
from sys import float_info
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Literal,
    Optional,
    TypedDict,
    TypeVar,
    Union,
)

from hypothesis.errors import HypothesisWarning
from hypothesis.internal.cache import LRUCache
from hypothesis.internal.compat import WINDOWS, int_from_bytes
from hypothesis.internal.conjecture.choice import (
    StringKWargs,
    choice_kwargs_key,
    choice_permitted,
)
from hypothesis.internal.conjecture.floats import float_to_lex, lex_to_float
from hypothesis.internal.conjecture.junkdrawer import bits_to_bytes
from hypothesis.internal.conjecture.utils import (
    INT_SIZES,
    INT_SIZES_SAMPLER,
    Sampler,
    many,
)
from hypothesis.internal.floats import (
    SIGNALING_NAN,
    float_to_int,
    make_float_clamper,
    next_down,
    next_up,
    sign_aware_lte,
)
from hypothesis.internal.intervalsets import IntervalSet

if TYPE_CHECKING:
    from typing import TypeAlias

    from hypothesis.internal.conjecture.data import ConjectureData

T = TypeVar("T")
_Lifetime: "TypeAlias" = Literal["test_case", "test_function"]
COLLECTION_DEFAULT_MAX_SIZE = 10**10  # "arbitrarily large"


# The available `PrimitiveProvider`s, and therefore also the available backends
# for use by @settings(backend=...). The key is the name to be used in the backend=
# value, and the value is the importable path to a subclass of PrimitiveProvider.
#
# See also
# https://hypothesis.readthedocs.io/en/latest/strategies.html#alternative-backends-for-hypothesis.
#
# NOTE: the PrimitiveProvider interface is not yet stable. We may continue to
# make breaking changes to it. (but if you want to experiment and don't mind
# breakage, here you go!)
AVAILABLE_PROVIDERS = {
    "hypothesis": "hypothesis.internal.conjecture.providers.HypothesisProvider",
    "hypothesis-urandom": "hypothesis.internal.conjecture.providers.URandomProvider",
}
FLOAT_INIT_LOGIC_CACHE = LRUCache(4096)
STRING_SAMPLER_CACHE = LRUCache(64)

NASTY_FLOATS = sorted(
    [
        0.0,
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
        9007199254740992,
        1 - 10e-6,
        2 + 10e-6,
        1.192092896e-07,
        2.2204460492503131e-016,
    ]
    + [2.0**-n for n in (24, 14, 149, 126)]  # minimum (sub)normals for float16,32
    + [float_info.min / n for n in (2, 10, 1000, 100_000)]  # subnormal in float64
    + [math.inf, math.nan] * 5
    + [SIGNALING_NAN],
    key=float_to_lex,
)
NASTY_FLOATS = list(map(float, NASTY_FLOATS))
NASTY_FLOATS.extend([-x for x in NASTY_FLOATS])

NASTY_STRINGS = sorted(
    [
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
        "æœÆŒﬀʤʨß"
        # emoticons
        "(╯°□°）╯︵ ┻━┻)",
        # emojis
        "😍",
        "🇺🇸",
        # emoji modifiers
        "🏻"  # U+1F3FB Light Skin Tone,
        "👍🏻",  # 👍 followed by U+1F3FB
        # RTL text
        "الكل في المجمو عة",
        # Ogham text, which contains the only character in the Space Separators
        # unicode category (Zs) that isn't visually blank:  .  # noqa: RUF003
        "᚛ᚄᚓᚐᚋᚒᚄ ᚑᚄᚂᚑᚏᚅ᚜",
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
    ],
    key=len,
)

# Masks for masking off the first byte of an n-bit buffer.
# The appropriate mask is stored at position n % 8.
BYTE_MASKS = [(1 << n) - 1 for n in range(8)]
BYTE_MASKS[0] = 255


class _BackendInfoMsg(TypedDict):
    type: str
    title: str
    content: Union[str, dict[str, Any]]


class PrimitiveProvider(abc.ABC):
    # This is the low-level interface which would also be implemented
    # by e.g. CrossHair, by an Atheris-hypothesis integration, etc.
    # We'd then build the structured tree handling, database and replay
    # support, etc. on top of this - so all backends get those for free.
    #
    # See https://github.com/HypothesisWorks/hypothesis/issues/3086

    # How long a provider instance is used for. One of test_function or
    # test_case. Defaults to test_function.
    #
    # If test_function, a single provider instance will be instantiated and used
    # for the entirety of each test function. I.e., roughly one provider per
    # @given annotation. This can be useful if you need to track state over many
    # executions to a test function.
    #
    # This lifetime will cause None to be passed for the ConjectureData object
    # in PrimitiveProvider.__init__, because that object is instantiated per
    # test case.
    #
    # If test_case, a new provider instance will be instantiated and used each
    # time hypothesis tries to generate a new input to the test function. This
    # lifetime can access the passed ConjectureData object.
    #
    # Non-hypothesis providers probably want to set a lifetime of test_function.
    lifetime: _Lifetime = "test_function"

    # Solver-based backends such as hypothesis-crosshair use symbolic values
    # which record operations performed on them in order to discover new paths.
    # If avoid_realization is set to True, hypothesis will avoid interacting with
    # symbolic choices returned by the provider in any way that would force the
    # solver to narrow the range of possible values for that symbolic.
    #
    # Setting this to True disables some hypothesis features, such as
    # DataTree-based deduplication, and some internal optimizations, such as
    # caching kwargs. Only enable this if it is necessary for your backend.
    avoid_realization = False

    def __init__(self, conjecturedata: Optional["ConjectureData"], /) -> None:
        self._cd = conjecturedata

    def per_test_case_context_manager(self):
        return contextlib.nullcontext()

    def realize(self, value: T) -> T:
        """
        Called whenever hypothesis requires a concrete (non-symbolic) value from
        a potentially symbolic value. Hypothesis will not check that `value` is
        symbolic before calling `realize`, so you should handle the case where
        `value` is non-symbolic.

        The returned value should be non-symbolic.  If you cannot provide a value,
        raise hypothesis.errors.BackendCannotProceed("discard_test_case")
        """
        return value

    def observe_test_case(self) -> dict[str, Any]:
        """Called at the end of the test case when observability mode is active.

        The return value should be a non-symbolic json-encodable dictionary,
        and will be included as `observation["metadata"]["backend"]`.
        """
        return {}

    def observe_information_messages(
        self, *, lifetime: _Lifetime
    ) -> Iterable[_BackendInfoMsg]:
        """Called at the end of each test case and again at end of the test function.

        Return an iterable of `{type: info/alert/error, title: str, content: str|dict}`
        dictionaries to be delivered as individual information messages.
        (Hypothesis adds the `run_start` timestamp and `property` name for you.)
        """
        assert lifetime in ("test_case", "test_function")
        yield from []

    @abc.abstractmethod
    def draw_boolean(
        self,
        p: float = 0.5,
    ) -> bool:
        raise NotImplementedError

    @abc.abstractmethod
    def draw_integer(
        self,
        min_value: Optional[int] = None,
        max_value: Optional[int] = None,
        *,
        # weights are for choosing an element index from a bounded range
        weights: Optional[dict[int, float]] = None,
        shrink_towards: int = 0,
    ) -> int:
        raise NotImplementedError

    @abc.abstractmethod
    def draw_float(
        self,
        *,
        min_value: float = -math.inf,
        max_value: float = math.inf,
        allow_nan: bool = True,
        smallest_nonzero_magnitude: float,
        # TODO: consider supporting these float widths at the IR level in the
        # future.
        # width: Literal[16, 32, 64] = 64,
        # exclude_min and exclude_max handled higher up,
    ) -> float:
        raise NotImplementedError

    @abc.abstractmethod
    def draw_string(
        self,
        intervals: IntervalSet,
        *,
        min_size: int = 0,
        max_size: int = COLLECTION_DEFAULT_MAX_SIZE,
    ) -> str:
        raise NotImplementedError

    @abc.abstractmethod
    def draw_bytes(
        self,
        min_size: int = 0,
        max_size: int = COLLECTION_DEFAULT_MAX_SIZE,
    ) -> bytes:
        raise NotImplementedError

    def span_start(self, label: int, /) -> None:  # noqa: B027  # non-abstract noop
        """Marks the beginning of a semantically meaningful span.

        Providers can optionally track this data to learn which sub-sequences
        of draws correspond to a higher-level object, recovering the parse tree.
        `label` is an opaque integer, which will be shared by all spans drawn
        from a particular strategy.

        This method is called from ConjectureData.start_span().
        """

    def span_end(self, discard: bool, /) -> None:  # noqa: B027, FBT001
        """Marks the end of a semantically meaningful span.

        `discard` is True when the draw was filtered out or otherwise marked as
        unlikely to contribute to the input data as seen by the user's test.
        Note however that side effects can make this determination unsound.

        This method is called from ConjectureData.stop_span().
        """


class HypothesisProvider(PrimitiveProvider):
    lifetime = "test_case"

    def __init__(self, conjecturedata: Optional["ConjectureData"], /):
        super().__init__(conjecturedata)
        self._random = None if self._cd is None else self._cd._random

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
        min_value: Optional[int] = None,
        max_value: Optional[int] = None,
        *,
        weights: Optional[dict[int, float]] = None,
        shrink_towards: int = 0,
    ) -> int:
        assert self._cd is not None

        center = 0
        if min_value is not None:
            center = max(min_value, center)
        if max_value is not None:
            center = min(max_value, center)

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
                return self._draw_bounded_integer(min_value, max_value)
            # implicit reliance on dicts being sorted for determinism
            return list(weights)[idx - 1]

        if min_value is None and max_value is None:
            return self._draw_unbounded_integer()

        if min_value is None:
            assert max_value is not None
            probe = max_value + 1
            while max_value < probe:
                probe = center + self._draw_unbounded_integer()
            return probe

        if max_value is None:
            assert min_value is not None
            probe = min_value - 1
            while probe < min_value:
                probe = center + self._draw_unbounded_integer()
            return probe

        return self._draw_bounded_integer(min_value, max_value)

    def draw_float(
        self,
        *,
        min_value: float = -math.inf,
        max_value: float = math.inf,
        allow_nan: bool = True,
        smallest_nonzero_magnitude: float,
        # TODO: consider supporting these float widths at the IR level in the
        # future.
        # width: Literal[16, 32, 64] = 64,
        # exclude_min and exclude_max handled higher up,
    ) -> float:
        (
            sampler,
            clamper,
            nasty_floats,
        ) = self._draw_float_init_logic(
            min_value=min_value,
            max_value=max_value,
            allow_nan=allow_nan,
            smallest_nonzero_magnitude=smallest_nonzero_magnitude,
        )

        assert self._cd is not None

        while True:
            i = sampler.sample(self._cd) if sampler else 0
            if i == 0:
                result = self._draw_float()
                if allow_nan and math.isnan(result):
                    clamped = result  # pragma: no cover
                else:
                    clamped = clamper(result)
                if float_to_int(clamped) != float_to_int(result) and not (
                    math.isnan(result) and allow_nan
                ):
                    result = clamped
            else:
                result = nasty_floats[i - 1]
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

        sampler, nasty_strings = self._draw_string_sampler(
            intervals=intervals,
            min_size=min_size,
            max_size=max_size,
        )

        if sampler is not None and self.draw_boolean(p=0.05):
            return nasty_strings[sampler.sample(self._cd)]

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

    def _draw_unbounded_integer(self) -> int:
        assert self._cd is not None
        assert self._random is not None

        size = INT_SIZES[INT_SIZES_SAMPLER.sample(self._cd)]

        r = self._random.getrandbits(size)
        sign = r & 1
        r >>= 1
        if sign:
            r = -r
        return r

    def _draw_bounded_integer(
        self,
        lower: int,
        upper: int,
        *,
        vary_size: bool = True,
    ) -> int:
        assert lower <= upper
        assert self._cd is not None
        assert self._random is not None

        if lower == upper:
            return lower

        bits = (upper - lower).bit_length()
        if bits > 24 and vary_size and self._random.random() < 7 / 8:
            # For large ranges, we combine the uniform random distribution
            # with a weighting scheme with moderate chance.  Cutoff at 2 ** 24 so that our
            # choice of unicode characters is uniform but the 32bit distribution is not.
            idx = INT_SIZES_SAMPLER.sample(self._cd)
            cap_bits = min(bits, INT_SIZES[idx])
            upper = min(upper, lower + 2**cap_bits - 1)
            return self._random.randint(lower, upper)

        return self._random.randint(lower, upper)

    @classmethod
    def _draw_float_init_logic(
        cls,
        *,
        min_value: float,
        max_value: float,
        allow_nan: bool,
        smallest_nonzero_magnitude: float,
    ) -> tuple[
        Optional[Sampler],
        Callable[[float], float],
        list[float],
    ]:
        """
        Caches initialization logic for draw_float, as an alternative to
        computing this for *every* float draw.
        """
        # float_to_int allows us to distinguish between e.g. -0.0 and 0.0,
        # even in light of hash(-0.0) == hash(0.0) and -0.0 == 0.0.
        key = (
            float_to_int(min_value),
            float_to_int(max_value),
            allow_nan,
            float_to_int(smallest_nonzero_magnitude),
        )
        if key in FLOAT_INIT_LOGIC_CACHE:
            return FLOAT_INIT_LOGIC_CACHE[key]

        result = cls._compute_draw_float_init_logic(
            min_value=min_value,
            max_value=max_value,
            allow_nan=allow_nan,
            smallest_nonzero_magnitude=smallest_nonzero_magnitude,
        )
        FLOAT_INIT_LOGIC_CACHE[key] = result
        return result

    @staticmethod
    def _compute_draw_float_init_logic(
        *,
        min_value: float,
        max_value: float,
        allow_nan: bool,
        smallest_nonzero_magnitude: float,
    ) -> tuple[
        Optional[Sampler],
        Callable[[float], float],
        list[float],
    ]:
        if smallest_nonzero_magnitude == 0.0:  # pragma: no cover
            raise FloatingPointError(
                "Got allow_subnormal=True, but we can't represent subnormal floats "
                "right now, in violation of the IEEE-754 floating-point "
                "specification.  This is usually because something was compiled with "
                "-ffast-math or a similar option, which sets global processor state.  "
                "See https://simonbyrne.github.io/notes/fastmath/ for a more detailed "
                "writeup - and good luck!"
            )

        def permitted(f: float) -> bool:
            if math.isnan(f):
                return allow_nan
            if 0 < abs(f) < smallest_nonzero_magnitude:
                return False
            return sign_aware_lte(min_value, f) and sign_aware_lte(f, max_value)

        boundary_values = [
            min_value,
            next_up(min_value),
            min_value + 1,
            max_value - 1,
            next_down(max_value),
            max_value,
        ]
        nasty_floats = [f for f in NASTY_FLOATS + boundary_values if permitted(f)]
        weights = [0.2 * len(nasty_floats)] + [0.8] * len(nasty_floats)
        sampler = Sampler(weights, observe=False) if nasty_floats else None

        clamper = make_float_clamper(
            min_value,
            max_value,
            smallest_nonzero_magnitude=smallest_nonzero_magnitude,
            allow_nan=allow_nan,
        )
        return (sampler, clamper, nasty_floats)

    @classmethod
    def _draw_string_sampler(
        cls,
        *,
        intervals: IntervalSet,
        min_size: int,
        max_size: int,
    ) -> tuple[Optional[Sampler], list[str]]:
        kwargs: StringKWargs = {
            "intervals": intervals,
            "min_size": min_size,
            "max_size": max_size,
        }
        key = choice_kwargs_key("string", kwargs)
        if key in STRING_SAMPLER_CACHE:
            return STRING_SAMPLER_CACHE[key]

        nasty_strings = [s for s in NASTY_STRINGS if choice_permitted(s, kwargs)]
        sampler = (
            Sampler([1 / len(nasty_strings)] * len(nasty_strings), observe=False)
            if nasty_strings
            else None
        )
        result = (sampler, nasty_strings)
        STRING_SAMPLER_CACHE[key] = result
        return result


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
        min_value: Optional[int] = None,
        max_value: Optional[int] = None,
        *,
        weights: Optional[dict[int, float]] = None,
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
        with open("/dev/urandom", "rb") as f:
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
