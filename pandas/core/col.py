from __future__ import annotations

from collections.abc import (
    Callable,
    Hashable,
)
from typing import (
    TYPE_CHECKING,
    Any,
)

from pandas.core.series import Series

if TYPE_CHECKING:
    from pandas import DataFrame


def parse_args(df: DataFrame, *args: Any) -> tuple[Series]:
    return tuple([x._func(df) if isinstance(x, Expr) else x for x in args])


def parse_kwargs(df: DataFrame, **kwargs: Any) -> dict[Hashable, Series]:
    return {
        key: val._func(df) if isinstance(val, Expr) else val
        for key, val in kwargs.items()
    }


class Expr:
    def __init__(self, func: Callable[[DataFrame], Any]) -> None:
        self._func = func

    def __call__(self, df: DataFrame) -> Series:
        result = self._func(df)
        if not isinstance(result, Series):
            msg = (
                "Expected function which returns Series, "
                f"got function which returns: {type(result)}"
            )
            raise TypeError(msg)
        return result

    # namespaces
    @property
    def dt(self) -> NamespaceExpr:
        return NamespaceExpr(self, "dt")

    @property
    def str(self) -> NamespaceExpr:
        return NamespaceExpr(self, "str")

    @property
    def cat(self) -> NamespaceExpr:
        return NamespaceExpr(self, "cat")

    @property
    def list(self) -> NamespaceExpr:
        return NamespaceExpr(self, "list")

    @property
    def sparse(self) -> NamespaceExpr:
        return NamespaceExpr(self, "sparse")

    @property
    def struct(self) -> NamespaceExpr:
        return NamespaceExpr(self, "struct")

    # Binary ops

    def __add__(self, other: Any) -> Expr:
        if isinstance(other, Expr):
            return Expr(lambda df: self._func(df).__add__(other._func(df)))
        return Expr(lambda df: self._func(df).__add__(other))

    def __radd__(self, other: Any) -> Expr:
        if isinstance(other, Expr):
            return Expr(lambda df: self._func(df).__radd__(other._func(df)))
        return Expr(lambda df: self._func(df).__radd__(other))

    def __sub__(self, other: Any) -> Expr:
        if isinstance(other, Expr):
            return Expr(lambda df: self._func(df).__sub__(other._func(df)))
        return Expr(lambda df: self._func(df).__sub__(other))

    def __rsub__(self, other: Any) -> Expr:
        if isinstance(other, Expr):
            return Expr(lambda df: self._func(df).__rsub__(other._func(df)))
        return Expr(lambda df: self._func(df).__rsub__(other))

    def __mul__(self, other: Any) -> Expr:
        if isinstance(other, Expr):
            return Expr(lambda df: self._func(df).__mul__(other._func(df)))
        return Expr(lambda df: self._func(df).__mul__(other))

    def __rmul__(self, other: Any) -> Expr:
        if isinstance(other, Expr):
            return Expr(lambda df: self._func(df).__rmul__(other._func(df)))
        return Expr(lambda df: self._func(df).__rmul__(other))

    def __truediv__(self, other: Any) -> Expr:
        if isinstance(other, Expr):
            return Expr(lambda df: self._func(df).__truediv__(other._func(df)))
        return Expr(lambda df: self._func(df).__truediv__(other))

    def __rtruediv__(self, other: Any) -> Expr:
        if isinstance(other, Expr):
            return Expr(lambda df: self._func(df).__rtruediv__(other._func(df)))
        return Expr(lambda df: self._func(df).__rtruediv__(other))

    def __floordiv__(self, other: Any) -> Expr:
        if isinstance(other, Expr):
            return Expr(lambda df: self._func(df).__floordiv__(other._func(df)))
        return Expr(lambda df: self._func(df).__floordiv__(other))

    def __rfloordiv__(self, other: Any) -> Expr:
        if isinstance(other, Expr):
            return Expr(lambda df: self._func(df).__rfloordiv__(other._func(df)))
        return Expr(lambda df: self._func(df).__rfloordiv__(other))

    def __ge__(self, other: Any) -> Expr:
        if isinstance(other, Expr):
            return Expr(lambda df: self._func(df).__ge__(other._func(df)))
        return Expr(lambda df: self._func(df).__ge__(other))

    def __gt__(self, other: Any) -> Expr:
        if isinstance(other, Expr):
            return Expr(lambda df: self._func(df).__gt__(other._func(df)))
        return Expr(lambda df: self._func(df).__gt__(other))

    def __le__(self, other: Any) -> Expr:
        if isinstance(other, Expr):
            return Expr(lambda df: self._func(df).__le__(other._func(df)))
        return Expr(lambda df: self._func(df).__le__(other))

    def __lt__(self, other: Any) -> Expr:
        if isinstance(other, Expr):
            return Expr(lambda df: self._func(df).__lt__(other._func(df)))
        return Expr(lambda df: self._func(df).__lt__(other))

    def __eq__(self, other: object) -> Expr:  # type: ignore[override]
        if isinstance(other, Expr):
            return Expr(lambda df: self._func(df).__eq__(other._func(df)))
        return Expr(lambda df: self._func(df).__eq__(other))

    def __ne__(self, other: object) -> Expr:  # type: ignore[override]
        if isinstance(other, Expr):
            return Expr(lambda df: self._func(df).__ne__(other._func(df)))
        return Expr(lambda df: self._func(df).__ne__(other))

    def __mod__(self, other: Any) -> Expr:
        if isinstance(other, Expr):
            return Expr(lambda df: self._func(df).__mod__(other._func(df)))
        return Expr(lambda df: self._func(df).__mod__(other))

    # Everything else

    # Function "pandas.core.col.Expr.str" is not valid as a type
    def __getattr__(self, attr: str, /) -> Any:  # type: ignore[valid-type]
        def func(df: DataFrame, *args: Any, **kwargs: Any) -> Any:
            parsed_args = parse_args(df, *args)
            parsed_kwargs = parse_kwargs(df, **kwargs)
            return getattr(self(df), attr)(*parsed_args, **parsed_kwargs)

        return lambda *args, **kwargs: Expr(lambda df: func(df, *args, **kwargs))


class NamespaceExpr:
    def __init__(self, func: Callable[[DataFrame], Any], namespace: str) -> None:
        self._func = func
        self._namespace = namespace

    def __getattr__(self, attr: str) -> Any:
        if isinstance(getattr(getattr(Series, self._namespace), attr), property):
            return Expr(
                lambda df: getattr(getattr(self._func(df), self._namespace), attr)
            )

        def func(df: DataFrame, *args: Any, **kwargs: Any) -> Any:
            parsed_args = parse_args(df, *args)
            parsed_kwargs = parse_kwargs(df, **kwargs)
            return getattr(getattr(self._func(df), self._namespace), attr)(
                *parsed_args, **parsed_kwargs
            )

        return lambda *args, **kwargs: Expr(lambda df: func(df, *args, **kwargs))


def col(col_name: Hashable) -> Expr:
    if not isinstance(col_name, Hashable):
        msg = f"Expected Hashable, got: {type(col_name)}"
        raise TypeError(msg)
    return Expr(lambda df: df[col_name])
