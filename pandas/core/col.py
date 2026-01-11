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

if TYPE_CHECKING:
    from pandas import (
        DataFrame,
        Series,
    )


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
    return tuple(
        [x._eval_expression(df) if isinstance(x, Expression) else x for x in args]
    )


def _parse_kwargs(df: DataFrame, **kwargs: Any) -> dict[str, Any]:
    # Parse `kwargs`, evaluating any expressions we encounter.
    return {
        key: val._eval_expression(df) if isinstance(val, Expression) else val
        for key, val in kwargs.items()
    }


def _pretty_print_args_kwargs(*args: Any, **kwargs: Any) -> str:
    inputs_repr = ", ".join(repr(arg) for arg in args)
    kwargs_repr = ", ".join(f"{k}={v!r}" for k, v in kwargs.items())

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

    def __init__(
        self,
        func: Callable[[DataFrame], Any],
        repr_str: str,
        needs_parenthese: bool = False,
    ) -> None:
        self._func = func
        self._repr_str = repr_str
        self._needs_parentheses = needs_parenthese

    def _eval_expression(self, df: DataFrame) -> Any:
        return self._func(df)

    def _with_op(
        self, op: str, other: Any, repr_str: str, needs_parentheses: bool = True
    ) -> Expression:
        if isinstance(other, Expression):
            return Expression(
                lambda df: getattr(self._eval_expression(df), op)(
                    other._eval_expression(df)
                ),
                repr_str,
                needs_parenthese=needs_parentheses,
            )
        else:
            return Expression(
                lambda df: getattr(self._eval_expression(df), op)(other),
                repr_str,
                needs_parenthese=needs_parentheses,
            )

    def _maybe_wrap_parentheses(self, other: Any) -> tuple[str, str]:
        if self._needs_parentheses:
            self_repr = f"({self!r})"
        else:
            self_repr = f"{self!r}"
        if isinstance(other, Expression) and other._needs_parentheses:
            other_repr = f"({other!r})"
        else:
            other_repr = f"{other!r}"
        return self_repr, other_repr

    # Binary ops
    def __add__(self, other: Any) -> Expression:
        self_repr, other_repr = self._maybe_wrap_parentheses(other)
        return self._with_op("__add__", other, f"{self_repr} + {other_repr}")

    def __radd__(self, other: Any) -> Expression:
        self_repr, other_repr = self._maybe_wrap_parentheses(other)
        return self._with_op("__radd__", other, f"{other_repr} + {self_repr}")

    def __sub__(self, other: Any) -> Expression:
        self_repr, other_repr = self._maybe_wrap_parentheses(other)
        return self._with_op("__sub__", other, f"{self_repr} - {other_repr}")

    def __rsub__(self, other: Any) -> Expression:
        self_repr, other_repr = self._maybe_wrap_parentheses(other)
        return self._with_op("__rsub__", other, f"{other_repr} - {self_repr}")

    def __mul__(self, other: Any) -> Expression:
        self_repr, other_repr = self._maybe_wrap_parentheses(other)
        return self._with_op("__mul__", other, f"{self_repr} * {other_repr}")

    def __rmul__(self, other: Any) -> Expression:
        self_repr, other_repr = self._maybe_wrap_parentheses(other)
        return self._with_op("__rmul__", other, f"{other_repr} * {self_repr}")

    def __truediv__(self, other: Any) -> Expression:
        self_repr, other_repr = self._maybe_wrap_parentheses(other)
        return self._with_op("__truediv__", other, f"{self_repr} / {other_repr}")

    def __rtruediv__(self, other: Any) -> Expression:
        self_repr, other_repr = self._maybe_wrap_parentheses(other)
        return self._with_op("__rtruediv__", other, f"{other_repr} / {self_repr}")

    def __floordiv__(self, other: Any) -> Expression:
        self_repr, other_repr = self._maybe_wrap_parentheses(other)
        return self._with_op("__floordiv__", other, f"{self_repr} // {other_repr}")

    def __rfloordiv__(self, other: Any) -> Expression:
        self_repr, other_repr = self._maybe_wrap_parentheses(other)
        return self._with_op("__rfloordiv__", other, f"{other_repr} // {self_repr}")

    def __ge__(self, other: Any) -> Expression:
        self_repr, other_repr = self._maybe_wrap_parentheses(other)
        return self._with_op("__ge__", other, f"{self_repr} >= {other_repr}")

    def __gt__(self, other: Any) -> Expression:
        self_repr, other_repr = self._maybe_wrap_parentheses(other)
        return self._with_op("__gt__", other, f"{self_repr} > {other_repr}")

    def __le__(self, other: Any) -> Expression:
        self_repr, other_repr = self._maybe_wrap_parentheses(other)
        return self._with_op("__le__", other, f"{self_repr} <= {other_repr}")

    def __lt__(self, other: Any) -> Expression:
        self_repr, other_repr = self._maybe_wrap_parentheses(other)
        return self._with_op("__lt__", other, f"{self_repr} < {other_repr}")

    def __eq__(self, other: object) -> Expression:  # type: ignore[override]
        self_repr, other_repr = self._maybe_wrap_parentheses(other)
        return self._with_op("__eq__", other, f"{self_repr} == {other_repr}")

    def __ne__(self, other: object) -> Expression:  # type: ignore[override]
        self_repr, other_repr = self._maybe_wrap_parentheses(other)
        return self._with_op("__ne__", other, f"{self_repr} != {other_repr}")

    def __mod__(self, other: Any) -> Expression:
        self_repr, other_repr = self._maybe_wrap_parentheses(other)
        return self._with_op("__mod__", other, f"{self_repr} % {other_repr}")

    def __rmod__(self, other: Any) -> Expression:
        self_repr, other_repr = self._maybe_wrap_parentheses(other)
        return self._with_op("__rmod__", other, f"{other_repr} % {self_repr}")

    # Logical ops
    def __and__(self, other: Any) -> Expression:
        self_repr, other_repr = self._maybe_wrap_parentheses(other)
        return self._with_op("__and__", other, f"{self_repr} & {other_repr}")

    def __rand__(self, other: Any) -> Expression:
        self_repr, other_repr = self._maybe_wrap_parentheses(other)
        return self._with_op("__rand__", other, f"{other_repr} & {self_repr}")

    def __or__(self, other: Any) -> Expression:
        self_repr, other_repr = self._maybe_wrap_parentheses(other)
        return self._with_op("__or__", other, f"{self_repr} | {other_repr}")

    def __ror__(self, other: Any) -> Expression:
        self_repr, other_repr = self._maybe_wrap_parentheses(other)
        return self._with_op("__ror__", other, f"{other_repr} | {self_repr}")

    def __xor__(self, other: Any) -> Expression:
        self_repr, other_repr = self._maybe_wrap_parentheses(other)
        return self._with_op("__xor__", other, f"{self_repr} ^ {other_repr}")

    def __rxor__(self, other: Any) -> Expression:
        self_repr, other_repr = self._maybe_wrap_parentheses(other)
        return self._with_op("__rxor__", other, f"{other_repr} ^ {self_repr}")

    def __invert__(self) -> Expression:
        return Expression(
            lambda df: ~self._eval_expression(df),
            f"~{self._repr_str}",
            needs_parenthese=True,
        )

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

    def __getitem__(self, item: Any) -> Expression:
        return self._with_op(
            "__getitem__", item, f"{self!r}[{item!r}]", needs_parentheses=True
        )

    def _call_with_func(self, func: Callable, **kwargs: Any) -> Expression:
        def wrapped(df: DataFrame) -> Any:
            parsed_kwargs = _parse_kwargs(df, **kwargs)
            return func(**parsed_kwargs)

        args_str = _pretty_print_args_kwargs(**kwargs)
        repr_str = func.__name__ + "(" + args_str + ")"

        return Expression(wrapped, repr_str)

    def __call__(self, *args: Any, **kwargs: Any) -> Expression:
        def func(df: DataFrame, *args: Any, **kwargs: Any) -> Any:
            parsed_args = _parse_args(df, *args)
            parsed_kwargs = _parse_kwargs(df, **kwargs)
            return self._eval_expression(df)(*parsed_args, **parsed_kwargs)

        args_str = _pretty_print_args_kwargs(*args, **kwargs)
        repr_str = f"{self._repr_str}({args_str})"
        return Expression(lambda df: func(df, *args, **kwargs), repr_str)

    def __getattr__(self, name: str, /) -> Any:
        repr_str = f"{self!r}"
        if self._needs_parentheses:
            repr_str = f"({repr_str})"
        repr_str += f".{name}"
        return Expression(lambda df: getattr(self._eval_expression(df), name), repr_str)

    def __repr__(self) -> str:
        return self._repr_str or "Expr(...)"


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
