from collections.abc import Callable, Generator, Iterable, Sequence
from types import ModuleType
from typing import Concatenate, Generic, override
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

###

type _IntegrandFunc[NumberT: npc.number] = Callable[Concatenate[onp.Array2D[NumberT], ...], onp.ArrayND[npc.number]]

_NumberT_co = TypeVar("_NumberT_co", bound=npc.number, default=np.float64, covariant=True)
_XPT_co = TypeVar("_XPT_co", default=ModuleType, covariant=True)

###

class Rule(Generic[_XPT_co]):  # undocumented
    xp: _XPT_co | None

    def estimate[NumberT: npc.number](
        self, /, f: _IntegrandFunc[NumberT], a: onp.ArrayND[NumberT], b: onp.ArrayND[NumberT], args: tuple[object, ...] = ()
    ) -> onp.ArrayND[NumberT]: ...  # abstract
    def estimate_error[NumberT: npc.number](
        self, /, f: _IntegrandFunc[NumberT], a: onp.ArrayND[NumberT], b: onp.ArrayND[NumberT], args: tuple[object, ...] = ()
    ) -> onp.ArrayND[NumberT]: ...

class FixedRule(Rule[_XPT_co], Generic[_XPT_co, _NumberT_co]):  # undocumented
    def __init__(self) -> None: ...
    @property
    def nodes_and_weights(self) -> tuple[onp.ArrayND[_NumberT_co], onp.ArrayND[_NumberT_co]]: ...  # abstract

class NestedFixedRule(FixedRule[_XPT_co, _NumberT_co], Generic[_XPT_co, _NumberT_co]):  # undocumented
    higher: FixedRule[_XPT_co, _NumberT_co]
    lower: FixedRule[_XPT_co, _NumberT_co]

    @override
    def __init__(self, /, higher: FixedRule[_XPT_co, _NumberT_co], lower: FixedRule[_XPT_co, _NumberT_co]) -> None: ...  # pyrefly:ignore[bad-override]
    @property
    def lower_nodes_and_weights(self) -> tuple[onp.ArrayND[_NumberT_co], onp.ArrayND[_NumberT_co]]: ...  # semi-abstract

class ProductNestedFixed(NestedFixedRule[_XPT_co, _NumberT_co], Generic[_XPT_co, _NumberT_co]):
    base_rules: Sequence[NestedFixedRule[_XPT_co, _NumberT_co]]
    @override
    def __init__(self, /, base_rules: Sequence[NestedFixedRule[_XPT_co, _NumberT_co]]) -> None: ...  # pyrefly:ignore[bad-override]

def _cartesian_product[NumberT: npc.number](arrays: Iterable[onp.ArrayND[NumberT]]) -> onp.Array2D[NumberT]: ...  # undocumented
def _split_subregion[NumberT: npc.number](
    a: onp.ArrayND[NumberT], b: onp.ArrayND[NumberT], xp: ModuleType, split_at: onp.ArrayND[NumberT] | None = None
) -> Generator[tuple[onp.ArrayND[NumberT], onp.ArrayND[NumberT]]]: ...  # undocumented
def _apply_fixed_rule[NumberT: npc.number](
    f: _IntegrandFunc[NumberT],
    a: onp.ArrayND[NumberT],
    b: onp.ArrayND[NumberT],
    orig_nodes: onp.ArrayND[npc.number],
    orig_weights: onp.ArrayND[npc.number],
    args: tuple[object, ...],
    xp: ModuleType,
) -> onp.ArrayND[NumberT]: ...  # undocumented
