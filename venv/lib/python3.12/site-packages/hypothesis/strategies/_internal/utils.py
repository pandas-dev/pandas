# This file is part of Hypothesis, which may be found at
# https://github.com/HypothesisWorks/hypothesis/
#
# Copyright the Hypothesis Authors.
# Individual contributors are listed in AUTHORS.rst and the git log.
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at https://mozilla.org/MPL/2.0/.

import dataclasses
import sys
from functools import partial
from inspect import signature
from typing import TYPE_CHECKING, Callable, TypeVar

import attr

from hypothesis.internal.cache import LRUReusedCache
from hypothesis.internal.floats import clamp, float_to_int
from hypothesis.internal.reflection import proxies
from hypothesis.vendor.pretty import pretty

if TYPE_CHECKING:
    from typing import TypeAlias

    from hypothesis.strategies._internal.strategies import SearchStrategy

T = TypeVar("T")
ValueKey: "TypeAlias" = tuple[type, object]
# (fn, args, kwargs)
StrategyCacheKey: "TypeAlias" = tuple[
    object, tuple[ValueKey, ...], frozenset[tuple[str, ValueKey]]
]

_strategies: dict[str, Callable[..., "SearchStrategy"]] = {}
# note: LRUReusedCache is already thread-local internally
_STRATEGY_CACHE = LRUReusedCache[StrategyCacheKey, object](1024)


def convert_value(v: object) -> ValueKey:
    if isinstance(v, float):
        return (float, float_to_int(v))
    return (type(v), v)


def clear_cache() -> None:
    _STRATEGY_CACHE.clear()


def cacheable(fn: T) -> T:
    from hypothesis.control import _current_build_context
    from hypothesis.strategies._internal.strategies import SearchStrategy

    @proxies(fn)
    def cached_strategy(*args, **kwargs):
        context = _current_build_context.value
        if context is not None and context.data.provider.avoid_realization:
            return fn(*args, **kwargs)

        try:
            kwargs_cache_key = {(k, convert_value(v)) for k, v in kwargs.items()}
        except TypeError:
            return fn(*args, **kwargs)
        cache_key = (fn, tuple(map(convert_value, args)), frozenset(kwargs_cache_key))
        try:
            if cache_key in _STRATEGY_CACHE:
                return _STRATEGY_CACHE[cache_key]
        except TypeError:
            return fn(*args, **kwargs)
        else:
            result = fn(*args, **kwargs)
            if not isinstance(result, SearchStrategy) or result.is_cacheable:
                result._is_singleton = True
                _STRATEGY_CACHE[cache_key] = result
            return result

    cached_strategy.__clear_cache = clear_cache  # type: ignore
    return cached_strategy


def defines_strategy(
    *,
    force_reusable_values: bool = False,
    try_non_lazy: bool = False,
    never_lazy: bool = False,
) -> Callable[[T], T]:
    """Returns a decorator for strategy functions.

    If ``force_reusable_values`` is True, the returned strategy will be marked
    with ``.has_reusable_values == True`` even if it uses maps/filters or
    non-reusable strategies internally. This tells our numpy/pandas strategies
    that they can implicitly use such strategies as background values.

    If ``try_non_lazy`` is True, attempt to execute the strategy definition
    function immediately, so that a LazyStrategy is only returned if this
    raises an exception.

    If ``never_lazy`` is True, the decorator performs no lazy-wrapping at all,
    and instead returns the original function.
    """

    def decorator(strategy_definition):
        """A decorator that registers the function as a strategy and makes it
        lazily evaluated."""
        _strategies[strategy_definition.__name__] = signature(strategy_definition)

        if never_lazy:
            assert not try_non_lazy
            # We could potentially support never_lazy + force_reusable_values
            # with a suitable wrapper, but currently there are no callers that
            # request this combination.
            assert not force_reusable_values
            return strategy_definition

        from hypothesis.strategies._internal.lazy import LazyStrategy

        @proxies(strategy_definition)
        def accept(*args, **kwargs):
            if try_non_lazy:
                # Why not try this unconditionally?  Because we'd end up with very
                # deep nesting of recursive strategies - better to be lazy unless we
                # *know* that eager evaluation is the right choice.
                try:
                    return strategy_definition(*args, **kwargs)
                except Exception:
                    # If invoking the strategy definition raises an exception,
                    # wrap that up in a LazyStrategy so it happens again later.
                    pass
            result = LazyStrategy(strategy_definition, args, kwargs)
            if force_reusable_values:
                # Setting `force_has_reusable_values` here causes the recursive
                # property code to set `.has_reusable_values == True`.
                result.force_has_reusable_values = True
                assert result.has_reusable_values
            return result

        accept.is_hypothesis_strategy_function = True
        return accept

    return decorator


def to_jsonable(obj: object, *, avoid_realization: bool) -> object:
    """Recursively convert an object to json-encodable form.

    This is not intended to round-trip, but rather provide an analysis-ready
    format for observability.  To avoid side affects, we pretty-print all but
    known types.
    """
    if isinstance(obj, (str, int, float, bool, type(None))):
        # We convert integers of 2**63 to floats, to avoid crashing external
        # utilities with a 64 bit integer cap (notable, sqlite). See
        # https://github.com/HypothesisWorks/hypothesis/pull/3797#discussion_r1413425110
        # and https://github.com/simonw/sqlite-utils/issues/605.
        if isinstance(obj, int) and not isinstance(obj, bool) and abs(obj) >= 2**63:
            # Silently clamp very large ints to max_float, to avoid OverflowError when
            # casting to float.  (but avoid adding more constraints to symbolic values)
            if avoid_realization:
                return "<symbolic>"
            obj = clamp(-sys.float_info.max, obj, sys.float_info.max)
            return float(obj)
        return obj
    if avoid_realization:
        return "<symbolic>"

    recur = partial(to_jsonable, avoid_realization=avoid_realization)
    if isinstance(obj, (list, tuple, set, frozenset)):
        if isinstance(obj, tuple) and hasattr(obj, "_asdict"):
            return recur(obj._asdict())  # treat namedtuples as dicts
        return [recur(x) for x in obj]
    if isinstance(obj, dict):
        return {
            k if isinstance(k, str) else pretty(k): recur(v) for k, v in obj.items()
        }

    # Hey, might as well try calling a .to_json() method - it works for Pandas!
    # We try this before the below general-purpose handlers to give folks a
    # chance to control this behavior on their custom classes.
    try:
        return recur(obj.to_json())  # type: ignore
    except Exception:
        pass

    # Special handling for dataclasses, attrs, and pydantic classes
    if (
        (dcs := sys.modules.get("dataclasses"))
        and dcs.is_dataclass(obj)
        and not isinstance(obj, type)
    ):
        # Avoid dataclasses.asdict here to ensure that inner to_json overrides
        # can get called as well
        return {
            field.name: recur(getattr(obj, field.name))
            for field in dataclasses.fields(obj)  # type: ignore
        }
    if attr.has(type(obj)):
        return recur(attr.asdict(obj, recurse=False))  # type: ignore
    if (pyd := sys.modules.get("pydantic")) and isinstance(obj, pyd.BaseModel):
        return recur(obj.model_dump())

    # If all else fails, we'll just pretty-print as a string.
    return pretty(obj)
