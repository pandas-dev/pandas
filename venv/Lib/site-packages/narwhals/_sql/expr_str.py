from __future__ import annotations

from typing import TYPE_CHECKING, Any, Generic

from narwhals._compliant import LazyExprNamespace
from narwhals._compliant.any_namespace import StringNamespace
from narwhals._sql.typing import SQLExprT
from narwhals._utils import is_pyspark_pre_4

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

    def contains(self, pattern: SQLExprT, *, literal: bool) -> SQLExprT:

        def func(expr: NativeExpr, pattern: NativeExpr) -> NativeExpr:
            func_name = "contains" if literal else "regexp_matches"
            return self._function(func_name, expr, pattern)

        compliant_pattern = (
            self.compliant.__narwhals_namespace__().lit(pattern, dtype=None)
            if isinstance(pattern, str)
            else pattern
        )
        return self.compliant._with_elementwise(func, pattern=compliant_pattern)

    def ends_with(self, suffix: SQLExprT) -> SQLExprT:
        return self.compliant._with_elementwise(
            lambda expr, suffix: self._function("ends_with", expr, suffix), suffix=suffix
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
                self._function("add", col_length, self._lit(offset + 1))
                if offset < 0
                else self._lit(offset + 1)
            )
            _length = self._lit(length) if length is not None else col_length
            return self._function("substr", expr, _offset, _length)

        return self.compliant._with_elementwise(func)

    def split(self, by: str) -> SQLExprT:
        # PySpark < 4.0's `split` expects a raw Python string for `pattern`,
        # not a Column literal.
        _is_pyspark_pre_4 = is_pyspark_pre_4(self.compliant._implementation)
        split_by = by if _is_pyspark_pre_4 else self._lit(by)
        return self.compliant._with_elementwise(
            lambda expr: self._function("str_split", expr, split_by)
        )

    def starts_with(self, prefix: SQLExprT) -> SQLExprT:
        return self.compliant._with_elementwise(
            lambda expr, prefix: self._function("starts_with", expr, prefix),
            prefix=prefix,
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

        # PySpark < 4.0's `lpad` expects raw Python values for `len` and `pad`,
        # not Column literals.
        _is_pyspark_pre_4 = is_pyspark_pre_4(self.compliant._implementation)

        def func(expr: NativeExpr) -> NativeExpr:
            less_than_width = self._function("length", expr) < self._lit(width)
            zero = "0" if _is_pyspark_pre_4 else self._lit("0")
            width_after_sign = width - 1 if _is_pyspark_pre_4 else self._lit(width - 1)
            full_width = width if _is_pyspark_pre_4 else self._lit(width)
            hyphen, plus = self._lit("-"), self._lit("+")

            starts_with_minus = self._function("starts_with", expr, hyphen)
            starts_with_plus = self._function("starts_with", expr, plus)
            substring = self._function("substr", expr, self._lit(2))
            padded_substring = self._function("lpad", substring, width_after_sign, zero)
            return self._when(
                self._function("and", starts_with_minus, less_than_width),
                self._function("concat", hyphen, padded_substring),
                self._when(
                    self._function("and", starts_with_plus, less_than_width),
                    self._function("concat", plus, padded_substring),
                    self._when(
                        less_than_width,
                        self._function("lpad", expr, full_width, zero),
                        expr,
                    ),
                ),
            )

        # can't use `_with_elementwise` due to `when` operator.
        # TODO(unassigned): implement `window_func` like we do in `Expr.cast`
        return self.compliant._with_callable(func)

    def pad_start(self, length: int, fill_char: str) -> SQLExprT:
        # PySpark < 4.0's `lpad` expects raw Python values for `len` and `pad`,
        # not Column literals.
        _is_pyspark_pre_4 = is_pyspark_pre_4(self.compliant._implementation)
        lpad_length = length if _is_pyspark_pre_4 else self._lit(length)
        lpad_fill = fill_char if _is_pyspark_pre_4 else self._lit(fill_char)

        def _pad_start(expr: NativeExpr) -> NativeExpr:
            return self._when(
                self._function("length", expr) < self._lit(length),
                self._function("lpad", expr, lpad_length, lpad_fill),
                expr,
            )

        return self.compliant._with_callable(_pad_start)

    def pad_end(self, length: int, fill_char: str) -> SQLExprT:
        # PySpark < 4.0's `rpad` expects raw Python values for `len` and `pad`,
        # not Column literals.
        _is_pyspark_pre_4 = is_pyspark_pre_4(self.compliant._implementation)
        rpad_length = length if _is_pyspark_pre_4 else self._lit(length)
        rpad_fill = fill_char if _is_pyspark_pre_4 else self._lit(fill_char)

        def _pad_end(expr: NativeExpr) -> NativeExpr:
            return self._when(
                self._function("length", expr) < self._lit(length),
                self._function("rpad", expr, rpad_length, rpad_fill),
                expr,
            )

        return self.compliant._with_callable(_pad_end)
