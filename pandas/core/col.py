from __future__ import annotations

from collections.abc import (
    Callable,
    Hashable,
)
from typing import (
    TYPE_CHECKING,
    Any,
)

from pandas.core.dtypes.common import is_scalar

from pandas.core.series import Series

if TYPE_CHECKING:
    from pandas import DataFrame


# Used only for generating the str repr of expressions.
OP_SYMBOLS = {
    "__add__": "+",
    "__radd__": "+",
    "__sub__": "-",
    "__rsub__": "-",
    "__mul__": "*",
    "__rmul__": "*",
    "__truediv__": "/",
    "__rtruediv__": "/",
    "__floordiv__": "//",
    "__rfloordiv__": "//",
    "__mod__": "%",
    "__rmod__": "%",
    "__ge__": ">=",
    "__gt__": ">",
    "__le__": "<=",
    "__lt__": "<",
    "__eq__": "==",
    "__ne__": "!=",
}


def parse_args(df: DataFrame, *args: Any) -> tuple[Series]:
    return tuple([x._func(df) if isinstance(x, Expr) else x for x in args])


def parse_kwargs(df: DataFrame, **kwargs: Any) -> dict[Hashable, Series]:
    return {
        key: val._func(df) if isinstance(val, Expr) else val
        for key, val in kwargs.items()
    }


class Expr:
    def __init__(
        self, func: Callable[[DataFrame], Any], repr_str: str | None = None
    ) -> None:
        self._func = func
        self._repr_str = repr_str

    def __call__(self, df: DataFrame) -> Series:
        result = self._func(df)
        if not (isinstance(result, Series) or is_scalar(result)):
            msg = (
                "Expected function which returns Series or scalar, "
                f"got function which returns: {type(result)}"
            )
            raise TypeError(msg)
        return result

    def _with_binary_op(self, op: str, other: Any) -> Expr:
        op_symbol = OP_SYMBOLS.get(op, op)

        if isinstance(other, Expr):
            if op.startswith("__r"):
                repr_str = f"({other._repr_str} {op_symbol} {self._repr_str})"
            else:
                repr_str = f"({self._repr_str} {op_symbol} {other._repr_str})"
            return Expr(
                lambda df: getattr(self._func(df), op)(other._func(df)), repr_str
            )
        else:
            if op.startswith("__r"):
                repr_str = f"({other!r} {op_symbol} {self._repr_str})"
            else:
                repr_str = f"({self._repr_str} {op_symbol} {other!r})"
            return Expr(lambda df: getattr(self._func(df), op)(other), repr_str)

    # Binary ops
    def __add__(self, other: Any) -> Expr:
        return self._with_binary_op("__add__", other)

    def __radd__(self, other: Any) -> Expr:
        return self._with_binary_op("__radd__", other)

    def __sub__(self, other: Any) -> Expr:
        return self._with_binary_op("__sub__", other)

    def __rsub__(self, other: Any) -> Expr:
        return self._with_binary_op("__rsub__", other)

    def __mul__(self, other: Any) -> Expr:
        return self._with_binary_op("__mul__", other)

    def __rmul__(self, other: Any) -> Expr:
        return self._with_binary_op("__rmul__", other)

    def __truediv__(self, other: Any) -> Expr:
        return self._with_binary_op("__truediv__", other)

    def __rtruediv__(self, other: Any) -> Expr:
        return self._with_binary_op("__rtruediv__", other)

    def __floordiv__(self, other: Any) -> Expr:
        return self._with_binary_op("__floordiv__", other)

    def __rfloordiv__(self, other: Any) -> Expr:
        return self._with_binary_op("__rfloordiv__", other)

    def __ge__(self, other: Any) -> Expr:
        return self._with_binary_op("__ge__", other)

    def __gt__(self, other: Any) -> Expr:
        return self._with_binary_op("__gt__", other)

    def __le__(self, other: Any) -> Expr:
        return self._with_binary_op("__le__", other)

    def __lt__(self, other: Any) -> Expr:
        return self._with_binary_op("__lt__", other)

    def __eq__(self, other: object) -> Expr:  # type: ignore[override]
        return self._with_binary_op("__eq__", other)

    def __ne__(self, other: object) -> Expr:  # type: ignore[override]
        return self._with_binary_op("__ne__", other)

    def __mod__(self, other: Any) -> Expr:
        return self._with_binary_op("__mod__", other)

    def __rmod__(self, other: Any) -> Expr:
        return self._with_binary_op("__rmod__", other)

    # Everything else
    def __getattr__(self, attr: str, /) -> Callable[..., Expr]:
        def func(df: DataFrame, *args: Any, **kwargs: Any) -> Any:
            parsed_args = parse_args(df, *args)
            parsed_kwargs = parse_kwargs(df, **kwargs)
            return getattr(self(df), attr)(*parsed_args, **parsed_kwargs)

        def wrapper(*args: Any, **kwargs: Any) -> Expr:
            # Create a readable representation for method calls
            args_repr = ", ".join(
                repr(arg._repr_str if isinstance(arg, Expr) else arg) for arg in args
            )
            kwargs_repr = ", ".join(
                f"{k}={v._repr_str if isinstance(v, Expr) else v!r}"
                for k, v in kwargs.items()
            )

            all_args = []
            if args_repr:
                all_args.append(args_repr)
            if kwargs_repr:
                all_args.append(kwargs_repr)

            args_str = ", ".join(all_args)
            repr_str = f"{self._repr_str}.{attr}({args_str})"

            return Expr(lambda df: func(df, *args, **kwargs), repr_str)

        return wrapper

    def __repr__(self) -> str:
        return self._repr_str or "Expr(...)"

    # Namespaces
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


class NamespaceExpr:
    def __init__(self, func: Expr, namespace: str) -> None:
        self._func = func
        self._namespace = namespace

    def __getattr__(self, attr: str) -> Any:
        if isinstance(getattr(getattr(Series, self._namespace), attr), property):
            repr_str = f"{self._func._repr_str}.{self._namespace}.{attr}"
            return Expr(
                lambda df: getattr(getattr(self._func(df), self._namespace), attr),
                repr_str,
            )

        def func(df: DataFrame, *args: Any, **kwargs: Any) -> Any:
            parsed_args = parse_args(df, *args)
            parsed_kwargs = parse_kwargs(df, **kwargs)
            return getattr(getattr(self._func(df), self._namespace), attr)(
                *parsed_args, **parsed_kwargs
            )

        def wrapper(*args: Any, **kwargs: Any) -> Expr:
            # Create a readable representation for namespace method calls
            args_repr = ", ".join(
                repr(arg._repr_str if isinstance(arg, Expr) else arg) for arg in args
            )
            kwargs_repr = ", ".join(
                f"{k}={v._repr_str if isinstance(v, Expr) else v!r}"
                for k, v in kwargs.items()
            )

            all_args = []
            if args_repr:
                all_args.append(args_repr)
            if kwargs_repr:
                all_args.append(kwargs_repr)

            args_str = ", ".join(all_args)
            repr_str = f"{self._func._repr_str}.{self._namespace}.{attr}({args_str})"

            return Expr(lambda df: func(df, *args, **kwargs), repr_str)

        return wrapper


def col(col_name: Hashable) -> Expr:
    if not isinstance(col_name, Hashable):
        msg = f"Expected Hashable, got: {type(col_name)}"
        raise TypeError(msg)

    def func(df: DataFrame) -> Series:
        if col_name not in df.columns:
            msg = (
                f"Column '{col_name}' not found in given DataFrame.\n\n"
                f"Hint: did you mean one of {df.columns.tolist()} instead?"
            )
            raise ValueError(msg)
        return df[col_name]

    return Expr(func, f"col({col_name!r})")
