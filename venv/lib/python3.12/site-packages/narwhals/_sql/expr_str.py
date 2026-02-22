from __future__ import annotations

import operator
from typing import TYPE_CHECKING, Any, Generic

from narwhals._compliant import LazyExprNamespace
from narwhals._compliant.any_namespace import StringNamespace
from narwhals._sql.typing import SQLExprT

if TYPE_CHECKING:
    from narwhals._compliant.expr import NativeExpr


class SQLExprStringNamespace(
    LazyExprNamespace[SQLExprT], StringNamespace[SQLExprT], Generic[SQLExprT]
):
    def _lit(self, value: Any) -> NativeExpr:
        return self.compliant._lit(value)  # type: ignore[no-any-return]

    def _function(self, name: str, *args: Any) -> NativeExpr:
        return self.compliant._function(name, *args)  # type: ignore[no-any-return]

    def _when(
        self, condition: Any, value: Any, otherwise: Any | None = None
    ) -> NativeExpr:
        return self.compliant._when(condition, value, otherwise)  # type: ignore[no-any-return]

    def contains(self, pattern: str, *, literal: bool) -> SQLExprT:
        def func(expr: NativeExpr) -> NativeExpr:
            if literal:
                return self._function("contains", expr, self._lit(pattern))
            return self._function("regexp_matches", expr, self._lit(pattern))

        return self.compliant._with_elementwise(func)

    def ends_with(self, suffix: str) -> SQLExprT:
        return self.compliant._with_elementwise(
            lambda expr: self._function("ends_with", expr, self._lit(suffix))
        )

    def len_chars(self) -> SQLExprT:
        return self.compliant._with_elementwise(
            lambda expr: self._function("length", expr)
        )

    def replace_all(self, value: SQLExprT, pattern: str, *, literal: bool) -> SQLExprT:
        fname: str = "replace" if literal else "regexp_replace"

        options: list[Any] = []
        if not literal and self.compliant._implementation.is_duckdb():
            options = [self._lit("g")]
        return self.compliant._with_elementwise(
            lambda expr, value: self._function(
                fname, expr, self._lit(pattern), value, *options
            ),
            value=value,
        )

    def slice(self, offset: int, length: int | None) -> SQLExprT:
        def func(expr: NativeExpr) -> NativeExpr:
            col_length = self._function("length", expr)

            _offset = (
                operator.add(col_length, self._lit(offset + 1))
                if offset < 0
                else self._lit(offset + 1)
            )
            _length = self._lit(length) if length is not None else col_length
            return self._function("substr", expr, _offset, _length)

        return self.compliant._with_elementwise(func)

    def split(self, by: str) -> SQLExprT:
        return self.compliant._with_elementwise(
            lambda expr: self._function("str_split", expr, self._lit(by))
        )

    def starts_with(self, prefix: str) -> SQLExprT:
        return self.compliant._with_elementwise(
            lambda expr: self._function("starts_with", expr, self._lit(prefix))
        )

    def strip_chars(self, characters: str | None) -> SQLExprT:
        import string

        return self.compliant._with_elementwise(
            lambda expr: self._function(
                "trim",
                expr,
                self._lit(string.whitespace if characters is None else characters),
            )
        )

    def to_lowercase(self) -> SQLExprT:
        return self.compliant._with_elementwise(
            lambda expr: self._function("lower", expr)
        )

    def to_uppercase(self) -> SQLExprT:
        return self.compliant._with_elementwise(
            lambda expr: self._function("upper", expr)
        )

    def zfill(self, width: int) -> SQLExprT:
        # There is no built-in zfill function, so we need to implement it manually
        # using string manipulation functions.

        def func(expr: NativeExpr) -> NativeExpr:
            less_than_width = self._function("length", expr) < self._lit(width)
            zero, hyphen, plus = self._lit("0"), self._lit("-"), self._lit("+")

            starts_with_minus = self._function("starts_with", expr, hyphen)
            starts_with_plus = self._function("starts_with", expr, plus)
            substring = self._function("substr", expr, self._lit(2))
            padded_substring = self._function(
                "lpad", substring, self._lit(width - 1), zero
            )
            return self._when(
                operator.and_(starts_with_minus, less_than_width),
                self._function("concat", hyphen, padded_substring),
                self._when(
                    operator.and_(starts_with_plus, less_than_width),
                    self._function("concat", plus, padded_substring),
                    self._when(
                        less_than_width,
                        self._function("lpad", expr, self._lit(width), zero),
                        expr,
                    ),
                ),
            )

        # can't use `_with_elementwise` due to `when` operator.
        # TODO(unassigned): implement `window_func` like we do in `Expr.cast`
        return self.compliant._with_callable(func)

    def pad_start(self, length: int, fill_char: str) -> SQLExprT:
        def _pad_start(expr: NativeExpr) -> NativeExpr:
            return self._when(
                self._function("length", expr) < self._lit(length),
                self._function("lpad", expr, self._lit(length), self._lit(fill_char)),
                expr,
            )

        return self.compliant._with_callable(_pad_start)

    def pad_end(self, length: int, fill_char: str) -> SQLExprT:
        def _pad_end(expr: NativeExpr) -> NativeExpr:
            return self._when(
                self._function("length", expr) < self._lit(length),
                self._function("rpad", expr, self._lit(length), self._lit(fill_char)),
                expr,
            )

        return self.compliant._with_callable(_pad_end)
