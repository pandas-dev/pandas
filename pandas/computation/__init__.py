import warnings


warnings.warn(
    "The pandas.computation module is deprecated and will be removed in a future "
    "version. Please import from the pandas.core.computation module instead.",
    FutureWarning, stacklevel=2
)


from . import (
    align, api, common, engines, eval,
    expr, expressions, ops, pytables, scope
)
