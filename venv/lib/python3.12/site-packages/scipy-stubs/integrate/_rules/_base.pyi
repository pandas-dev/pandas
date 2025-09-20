from collections.abc import Callable, Generator, Iterable, Sequence
from types import ModuleType
from typing import Concatenate, Generic, TypeAlias
from typing_extensions import TypeVar, override

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

_NumberT = TypeVar("_NumberT", bound=npc.number)
_NumberT_co = TypeVar("_NumberT_co", bound=npc.number, default=np.float64, covariant=True)
_XPT_co = TypeVar("_XPT_co", default=ModuleType, covariant=True)

_IntegrandFunc: TypeAlias = Callable[Concatenate[onp.Array2D[_NumberT], ...], onp.ArrayND[npc.number]]

###

class Rule(Generic[_XPT_co]):  # undocumented
    xp: _XPT_co | None

    def estimate(
        self, /, f: _IntegrandFunc[_NumberT], a: onp.ArrayND[_NumberT], b: onp.ArrayND[_NumberT], args: tuple[object, ...] = ()
    ) -> onp.ArrayND[_NumberT]: ...  # abstract
    def estimate_error(
        self, /, f: _IntegrandFunc[_NumberT], a: onp.ArrayND[_NumberT], b: onp.ArrayND[_NumberT], args: tuple[object, ...] = ()
    ) -> onp.ArrayND[_NumberT]: ...

class FixedRule(Rule[_XPT_co], Generic[_XPT_co, _NumberT_co]):  # undocumented
    def __init__(self) -> None: ...
    @property
    def nodes_and_weights(self) -> tuple[onp.ArrayND[_NumberT_co], onp.ArrayND[_NumberT_co]]: ...  # abstract

class NestedFixedRule(FixedRule[_XPT_co, _NumberT_co], Generic[_XPT_co, _NumberT_co]):  # undocumented
    higher: FixedRule[_XPT_co, _NumberT_co]
    lower: FixedRule[_XPT_co, _NumberT_co]
    @override
    def __init__(self, /, higher: FixedRule[_XPT_co, _NumberT_co], lower: FixedRule[_XPT_co, _NumberT_co]) -> None: ...
    @property
    def lower_nodes_and_weights(self) -> tuple[onp.ArrayND[_NumberT_co], onp.ArrayND[_NumberT_co]]: ...  # semi-abstract

class ProductNestedFixed(NestedFixedRule[_XPT_co, _NumberT_co], Generic[_XPT_co, _NumberT_co]):
    base_rules: Sequence[NestedFixedRule[_XPT_co, _NumberT_co]]
    @override
    def __init__(self, /, base_rules: Sequence[NestedFixedRule[_XPT_co, _NumberT_co]]) -> None: ...

def _cartesian_product(arrays: Iterable[onp.ArrayND[_NumberT]]) -> onp.Array2D[_NumberT]: ...  # undocumented
def _split_subregion(
    a: onp.ArrayND[_NumberT], b: onp.ArrayND[_NumberT], xp: ModuleType, split_at: onp.ArrayND[_NumberT] | None = None
) -> Generator[tuple[onp.ArrayND[_NumberT], onp.ArrayND[_NumberT]]]: ...  # undocumented
def _apply_fixed_rule(
    f: _IntegrandFunc[_NumberT],
    a: onp.ArrayND[_NumberT],
    b: onp.ArrayND[_NumberT],
    orig_nodes: onp.ArrayND[npc.number],
    orig_weights: onp.ArrayND[npc.number],
    args: tuple[object, ...],
    xp: ModuleType,
) -> onp.ArrayND[_NumberT]: ...  # undocumented
