# This file is part of Hypothesis, which may be found at
# https://github.com/HypothesisWorks/hypothesis/
#
# Copyright the Hypothesis Authors.
# Individual contributors are listed in AUTHORS.rst and the git log.
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at https://mozilla.org/MPL/2.0/.

from collections.abc import Hashable
from typing import Any, Optional

from hypothesis.internal.conjecture.data import ConjectureData
from hypothesis.strategies._internal import SearchStrategy
from hypothesis.strategies._internal.strategies import Ex


class SharedStrategy(SearchStrategy[Ex]):
    def __init__(self, base: SearchStrategy[Ex], key: Optional[Hashable] = None):
        self.key = key
        self.base = base

    @property
    def supports_find(self) -> bool:
        return self.base.supports_find

    def __repr__(self) -> str:
        if self.key is not None:
            return f"shared({self.base!r}, key={self.key!r})"
        else:
            return f"shared({self.base!r})"

    # Ideally would be -> Ex, but key collisions with different-typed values are
    # possible. See https://github.com/HypothesisWorks/hypothesis/issues/4301.
    def do_draw(self, data: ConjectureData) -> Any:
        if self.key is None or getattr(self.base, "_is_singleton", False):
            strat_label = id(self.base)
        else:
            # Assume that uncached strategies are distinguishable by their
            # label. False negatives (even collisions w/id above) are ok as
            # long as they are infrequent.
            strat_label = self.base.label
        key = self.key or self
        if key not in data._shared_strategy_draws:
            drawn = data.draw(self.base)
            data._shared_strategy_draws[key] = (strat_label, drawn)
        else:
            drawn_strat_label, drawn = data._shared_strategy_draws[key]
            # Check disabled pending resolution of #4301
            if drawn_strat_label != strat_label:  # pragma: no cover
                pass
                # warnings.warn(
                #     f"Different strategies are shared under {key=}. This"
                #     " risks drawing values that are not valid examples for the strategy,"
                #     " or that have a narrower range than expected.",
                #     HypothesisWarning,
                #     stacklevel=1,
                # )
        return drawn
