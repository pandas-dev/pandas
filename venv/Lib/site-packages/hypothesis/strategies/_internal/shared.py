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
        key = self.key or self
        if key not in data._shared_strategy_draws:
            data._shared_strategy_draws[key] = data.draw(self.base)
        return data._shared_strategy_draws[key]
