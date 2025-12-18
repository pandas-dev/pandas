from __future__ import annotations

from collections.abc import (
    Callable,
    Hashable,
)
from typing import (
    TYPE_CHECKING,
    Any,
)

from pandas.util._decorators import set_module

from pandas.core.series import Series

if TYPE_CHECKING:
    from pandas import DataFrame


# Used only for generating the str repr of expressions.
_OP_SYMBOLS = {
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
    "__and__": "&",
    "__rand__": "&",
    "__or__": "|",
    "__ror__": "|",
    "__xor__": "^",
    "__rxor__": "^",
}


def _parse_args(df: DataFrame, *args: Any) -> tuple[Series]:
    # Parse `args`, evaluating any expressions we encounter.
    return tuple([x(df) if isinstance(x, Expression) else x for x in args])


def _parse_kwargs(df: DataFrame, **kwargs: Any) -> dict[str, Any]:
    # Parse `kwargs`, evaluating any expressions we encounter.
    return {
        key: val(df) if isinstance(val, Expression) else val
        for key, val in kwargs.items()
    }


def _pretty_print_args_kwargs(*args: Any, **kwargs: Any) -> str:
    inputs_repr = ", ".join(
        arg._repr_str if isinstance(arg, Expression) else repr(arg) for arg in args
    )
    kwargs_repr = ", ".join(
        f"{k}={v._repr_str if isinstance(v, Expression) else v!r}"
        for k, v in kwargs.items()
    )

    all_args = []
    if inputs_repr:
        all_args.append(inputs_repr)
    if kwargs_repr:
        all_args.append(kwargs_repr)

    return ", ".join(all_args)


@set_module("pandas.api.typing")
class Expression:
    """
    Class representing a deferred column.

    This is not meant to be instantiated directly. Instead, use :meth:`pandas.col`.
    """

    def __init__(self, func: Callable[[DataFrame], Any], repr_str: str) -> None:
        self._func = func
        self._repr_str = repr_str

    def __call__(self, df: DataFrame) -> Any:
        return self._func(df)

    def _with_binary_op(self, op: str, other: Any) -> Expression:
        op_symbol = _OP_SYMBOLS.get(op, op)

        if isinstance(other, Expression):
            if op.startswith("__r"):
                repr_str = f"({other._repr_str} {op_symbol} {self._repr_str})"
            else:
                repr_str = f"({self._repr_str} {op_symbol} {other._repr_str})"
            return Expression(lambda df: getattr(self(df), op)(other(df)), repr_str)
        else:
            if op.startswith("__r"):
                repr_str = f"({other!r} {op_symbol} {self._repr_str})"
            else:
                repr_str = f"({self._repr_str} {op_symbol} {other!r})"
            return Expression(lambda df: getattr(self(df), op)(other), repr_str)

    # Binary ops
    def __add__(self, other: Any) -> Expression:
        return self._with_binary_op("__add__", other)

    def __radd__(self, other: Any) -> Expression:
        return self._with_binary_op("__radd__", other)

    def __sub__(self, other: Any) -> Expression:
        return self._with_binary_op("__sub__", other)

    def __rsub__(self, other: Any) -> Expression:
        return self._with_binary_op("__rsub__", other)

    def __mul__(self, other: Any) -> Expression:
        return self._with_binary_op("__mul__", other)

    def __rmul__(self, other: Any) -> Expression:
        return self._with_binary_op("__rmul__", other)

    def __truediv__(self, other: Any) -> Expression:
        return self._with_binary_op("__truediv__", other)

    def __rtruediv__(self, other: Any) -> Expression:
        return self._with_binary_op("__rtruediv__", other)

    def __floordiv__(self, other: Any) -> Expression:
        return self._with_binary_op("__floordiv__", other)

    def __rfloordiv__(self, other: Any) -> Expression:
        return self._with_binary_op("__rfloordiv__", other)

    def __ge__(self, other: Any) -> Expression:
        return self._with_binary_op("__ge__", other)

    def __gt__(self, other: Any) -> Expression:
        return self._with_binary_op("__gt__", other)

    def __le__(self, other: Any) -> Expression:
        return self._with_binary_op("__le__", other)

    def __lt__(self, other: Any) -> Expression:
        return self._with_binary_op("__lt__", other)

    def __eq__(self, other: object) -> Expression:  # type: ignore[override]
        return self._with_binary_op("__eq__", other)

    def __ne__(self, other: object) -> Expression:  # type: ignore[override]
        return self._with_binary_op("__ne__", other)

    def __mod__(self, other: Any) -> Expression:
        return self._with_binary_op("__mod__", other)

    def __rmod__(self, other: Any) -> Expression:
        return self._with_binary_op("__rmod__", other)

    # Logical ops
    def __and__(self, other: Any) -> Expression:
        return self._with_binary_op("__and__", other)

    def __rand__(self, other: Any) -> Expression:
        return self._with_binary_op("__rand__", other)

    def __or__(self, other: Any) -> Expression:
        return self._with_binary_op("__or__", other)

    def __ror__(self, other: Any) -> Expression:
        return self._with_binary_op("__ror__", other)

    def __xor__(self, other: Any) -> Expression:
        return self._with_binary_op("__xor__", other)

    def __rxor__(self, other: Any) -> Expression:
        return self._with_binary_op("__rxor__", other)

    def __invert__(self) -> Expression:
        return Expression(lambda df: ~self(df), f"(~{self._repr_str})")

    def __array_ufunc__(
        self, ufunc: Callable[..., Any], method: str, *inputs: Any, **kwargs: Any
    ) -> Expression:
        def func(df: DataFrame) -> Any:
            parsed_inputs = _parse_args(df, *inputs)
            parsed_kwargs = _parse_kwargs(df, *kwargs)
            return ufunc(*parsed_inputs, **parsed_kwargs)

        args_str = _pretty_print_args_kwargs(*inputs, **kwargs)
        repr_str = f"{ufunc.__name__}({args_str})"

        return Expression(func, repr_str)

    # Everything else
    def __getattr__(self, attr: str, /) -> Any:
        if attr in Series._accessors:
            return NamespaceExpression(self, attr)

        def func(df: DataFrame, *args: Any, **kwargs: Any) -> Any:
            parsed_args = _parse_args(df, *args)
            parsed_kwargs = _parse_kwargs(df, **kwargs)
            return getattr(self(df), attr)(*parsed_args, **parsed_kwargs)

        def wrapper(*args: Any, **kwargs: Any) -> Expression:
            args_str = _pretty_print_args_kwargs(*args, **kwargs)
            repr_str = f"{self._repr_str}.{attr}({args_str})"

            return Expression(lambda df: func(df, *args, **kwargs), repr_str)

        return wrapper

    def __repr__(self) -> str:
        return self._repr_str or "Expr(...)"


class NamespaceExpression:
    def __init__(self, func: Expression, namespace: str) -> None:
        self._func = func
        self._namespace = namespace

    def __call__(self, df: DataFrame) -> Any:
        return self._func(df)

    def __getattr__(self, attr: str) -> Any:
        if isinstance(getattr(getattr(Series, self._namespace), attr), property):
            repr_str = f"{self._func._repr_str}.{self._namespace}.{attr}"
            return Expression(
                lambda df: getattr(getattr(self(df), self._namespace), attr),
                repr_str,
            )

        def func(df: DataFrame, *args: Any, **kwargs: Any) -> Any:
            parsed_args = _parse_args(df, *args)
            parsed_kwargs = _parse_kwargs(df, **kwargs)
            return getattr(getattr(self(df), self._namespace), attr)(
                *parsed_args, **parsed_kwargs
            )

        def wrapper(*args: Any, **kwargs: Any) -> Expression:
            args_str = _pretty_print_args_kwargs(*args, **kwargs)
            repr_str = f"{self._func._repr_str}.{self._namespace}.{attr}({args_str})"
            return Expression(lambda df: func(df, *args, **kwargs), repr_str)

        return wrapper


@set_module("pandas")
def col(col_name: Hashable) -> Expression:
    """
    Generate deferred object representing a column of a DataFrame.

    Any place which accepts ``lambda df: df[col_name]``, such as
    :meth:`DataFrame.assign` or :meth:`DataFrame.loc`, can also accept
    ``pd.col(col_name)``.

    .. versionadded:: 3.0.0

    Parameters
    ----------
    col_name : Hashable
        Column name.

    Returns
    -------
    `pandas.api.typing.Expression`
        A deferred object representing a column of a DataFrame.

    See Also
    --------
    DataFrame.query : Query columns of a dataframe using string expressions.

    Examples
    --------

    You can use `col` in `assign`.

    >>> df = pd.DataFrame({"name": ["beluga", "narwhal"], "speed": [100, 110]})
    >>> df.assign(name_titlecase=pd.col("name").str.title())
          name  speed name_titlecase
    0   beluga    100         Beluga
    1  narwhal    110        Narwhal

    You can also use it for filtering.

    >>> df.loc[pd.col("speed") > 105]
          name  speed
    1  narwhal    110
    """
    if not isinstance(col_name, Hashable):
        msg = f"Expected Hashable, got: {type(col_name)}"
        raise TypeError(msg)

    def func(df: DataFrame) -> Series:
        if col_name not in df.columns:
            columns_str = str(df.columns.tolist())
            max_len = 90
            if len(columns_str) > max_len:
                columns_str = columns_str[:max_len] + "...]"

            msg = (
                f"Column '{col_name}' not found in given DataFrame.\n\n"
                f"Hint: did you mean one of {columns_str} instead?"
            )
            raise ValueError(msg)
        return df[col_name]

    return Expression(func, f"col({col_name!r})")


__all__ = ["Expression", "col"]
