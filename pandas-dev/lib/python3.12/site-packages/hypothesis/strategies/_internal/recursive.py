# This file is part of Hypothesis, which may be found at
# https://github.com/HypothesisWorks/hypothesis/
#
# Copyright the Hypothesis Authors.
# Individual contributors are listed in AUTHORS.rst and the git log.
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at https://mozilla.org/MPL/2.0/.

import threading
import warnings
from collections.abc import Callable
from contextlib import contextmanager

from hypothesis.errors import HypothesisWarning, InvalidArgument
from hypothesis.internal.reflection import (
    get_pretty_function_description,
    is_first_param_referenced_in_function,
    is_identity_function,
)
from hypothesis.internal.validation import check_type
from hypothesis.strategies._internal.strategies import (
    OneOfStrategy,
    SearchStrategy,
    check_strategy,
)
from hypothesis.utils.deprecation import note_deprecation


class LimitReached(BaseException):
    pass


class LimitedStrategy(SearchStrategy):
    def __init__(self, strategy):
        super().__init__()
        self.base_strategy = strategy
        self._threadlocal = threading.local()

    @property
    def marker(self):
        return getattr(self._threadlocal, "marker", 0)

    @marker.setter
    def marker(self, value):
        self._threadlocal.marker = value

    @property
    def currently_capped(self):
        return getattr(self._threadlocal, "currently_capped", False)

    @currently_capped.setter
    def currently_capped(self, value):
        self._threadlocal.currently_capped = value

    def __repr__(self) -> str:
        return f"LimitedStrategy({self.base_strategy!r})"

    def do_validate(self) -> None:
        self.base_strategy.validate()

    def do_draw(self, data):
        assert self.currently_capped
        if self.marker <= 0:
            raise LimitReached
        self.marker -= 1
        return data.draw(self.base_strategy)

    @contextmanager
    def capped(self, max_templates):
        try:
            was_capped = self.currently_capped
            self.currently_capped = True
            self.marker = max_templates
            yield
        finally:
            self.currently_capped = was_capped


class RecursiveStrategy(SearchStrategy):
    def __init__(
        self,
        base: SearchStrategy,
        extend: Callable[[SearchStrategy], SearchStrategy],
        min_leaves: int | None,
        max_leaves: int,
    ):
        super().__init__()
        self.min_leaves = min_leaves
        self.max_leaves = max_leaves
        self.base = base
        self.limited_base = LimitedStrategy(base)
        self.extend = extend

        strategies = [self.limited_base, self.extend(self.limited_base)]
        while 2 ** (len(strategies) - 1) <= max_leaves:
            strategies.append(extend(OneOfStrategy(tuple(strategies))))
        # If min_leaves > 1, we can never draw from base directly
        if min_leaves is not None and min_leaves > 1:
            strategies = strategies[1:]
        self.strategy = OneOfStrategy(strategies)

    def __repr__(self) -> str:
        if not hasattr(self, "_cached_repr"):
            self._cached_repr = (
                f"recursive({self.base!r}, "
                f"{get_pretty_function_description(self.extend)}, "
                f"min_leaves={self.min_leaves}, max_leaves={self.max_leaves})"
            )
        return self._cached_repr

    def do_validate(self) -> None:
        check_strategy(self.base, "base")
        extended = self.extend(self.limited_base)
        check_strategy(extended, f"extend({self.limited_base!r})")
        self.limited_base.validate()
        extended.validate()

        if is_identity_function(self.extend):
            warnings.warn(
                "extend=lambda x: x is a no-op; you probably want to use a "
                "different extend function, or just use the base strategy directly.",
                HypothesisWarning,
                stacklevel=5,
            )

        if not is_first_param_referenced_in_function(self.extend):
            msg = (
                f"extend={get_pretty_function_description(self.extend)} doesn't use "
                "it's argument, and thus can't actually recurse!"
            )
            if self.min_leaves is None:
                note_deprecation(
                    msg,
                    since="2026-01-12",
                    has_codemod=False,
                    stacklevel=1,
                )
            else:
                raise InvalidArgument(msg)

        if self.min_leaves is not None:
            check_type(int, self.min_leaves, "min_leaves")
        check_type(int, self.max_leaves, "max_leaves")
        if self.min_leaves is not None and self.min_leaves <= 0:
            raise InvalidArgument(
                f"min_leaves={self.min_leaves!r} must be greater than zero"
            )
        if self.max_leaves <= 0:
            raise InvalidArgument(
                f"max_leaves={self.max_leaves!r} must be greater than zero"
            )
        if (self.min_leaves or 1) > self.max_leaves:
            raise InvalidArgument(
                f"min_leaves={self.min_leaves!r} must be less than or equal to "
                f"max_leaves={self.max_leaves!r}"
            )

    def do_draw(self, data):
        min_leaves_retries = 0
        while True:
            try:
                with self.limited_base.capped(self.max_leaves):
                    result = data.draw(self.strategy)
                    leaves_drawn = self.max_leaves - self.limited_base.marker
                    if self.min_leaves and leaves_drawn < self.min_leaves:
                        data.events[
                            f"Draw for {self!r} had fewer than "
                            f"min_leaves={self.min_leaves} and had to be retried"
                        ] = ""
                        min_leaves_retries += 1
                        if min_leaves_retries < 5:
                            continue
                        data.mark_invalid(f"min_leaves={self.min_leaves} unsatisfied")
                    return result
            except LimitReached:
                data.events[
                    f"Draw for {self!r} exceeded "
                    f"max_leaves={self.max_leaves} and had to be retried"
                ] = ""
