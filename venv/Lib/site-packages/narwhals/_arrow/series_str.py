from __future__ import annotations

import string
from typing import TYPE_CHECKING

import pyarrow as pa
import pyarrow.compute as pc

from narwhals._arrow.utils import (
    ArrowSeriesNamespace,
    extract_native,
    lit,
    parse_datetime_format,
    parse_time_format,
)
from narwhals._compliant.any_namespace import StringNamespace

if TYPE_CHECKING:
    from narwhals._arrow.series import ArrowSeries
    from narwhals._arrow.typing import Incomplete


class ArrowSeriesStringNamespace(ArrowSeriesNamespace, StringNamespace["ArrowSeries"]):
    def len_chars(self) -> ArrowSeries:
        return self.with_native(pc.utf8_length(self.native))

    def replace(
        self, value: ArrowSeries, pattern: str, *, literal: bool, n: int
    ) -> ArrowSeries:
        fn = pc.replace_substring if literal else pc.replace_substring_regex
        _, value_native = extract_native(self.compliant, value)
        if not isinstance(value_native, pa.StringScalar):
            msg = "PyArrow backed `.str.replace` only supports str replacement values"
            raise TypeError(msg)
        arr = fn(
            self.native, pattern, replacement=value_native.as_py(), max_replacements=n
        )
        return self.with_native(arr)

    def replace_all(
        self, value: ArrowSeries, pattern: str, *, literal: bool
    ) -> ArrowSeries:
        return self.replace(value, pattern, literal=literal, n=-1)

    def strip_chars(self, characters: str | None) -> ArrowSeries:
        return self.with_native(
            pc.utf8_trim(self.native, characters or string.whitespace)
        )

    def starts_with(self, prefix: ArrowSeries) -> ArrowSeries:
        _, prefix_native = extract_native(self.compliant, prefix)
        if not isinstance(prefix_native, pa.StringScalar):
            msg = "`.str.starts_with` only supports str prefix values for pyarrow backend"
            raise TypeError(msg)
        return self.with_native(
            pc.equal(
                self.slice(0, len(prefix_native.as_py())).native,
                lit(prefix_native.as_py()),
            )
        )

    def ends_with(self, suffix: ArrowSeries) -> ArrowSeries:
        _, suffix_native = extract_native(self.compliant, suffix)
        if not isinstance(suffix_native, pa.StringScalar):
            msg = "`.str.ends_with` only supports str suffix values for pyarrow backend"
            raise TypeError(msg)
        return self.with_native(
            pc.equal(
                self.slice(-len(suffix_native.as_py()), None).native,
                lit(suffix_native.as_py()),
            )
        )

    def contains(self, pattern: ArrowSeries, *, literal: bool) -> ArrowSeries:
        _, pattern_native = extract_native(self.compliant, pattern)
        if not isinstance(pattern_native, pa.StringScalar):
            msg = "`.str.contains` only supports str pattern values for pyarrow backend"
            raise TypeError(msg)
        fn = pc.match_substring if literal else pc.match_substring_regex
        return self.with_native(fn(self.native, pattern_native.as_py()))

    def slice(self, offset: int, length: int | None) -> ArrowSeries:
        stop = offset + length if length is not None else None
        return self.with_native(
            pc.utf8_slice_codeunits(self.native, start=offset, stop=stop)
        )

    def split(self, by: str) -> ArrowSeries:
        split_series = pc.split_pattern(self.native, by)  # type: ignore[call-overload]
        return self.with_native(split_series)

    def to_datetime(self, format: str | None) -> ArrowSeries:
        format = parse_datetime_format(self.native) if format is None else format
        timestamp_array = pc.strptime(self.native, format=format, unit="us")
        return self.with_native(timestamp_array)

    def to_date(self, format: str | None) -> ArrowSeries:
        return self.to_datetime(format=format).dt.date()

    def to_time(self, format: str | None) -> ArrowSeries:
        format = parse_time_format(self.native) if format is None else format
        timestamp_array = pc.strptime(self.native, format=format, unit="us")

        nw_time_dtype = self.version.dtypes.Time()
        return self.with_native(timestamp_array).cast(nw_time_dtype)

    def to_uppercase(self) -> ArrowSeries:
        return self.with_native(pc.utf8_upper(self.native))

    def to_lowercase(self) -> ArrowSeries:
        return self.with_native(pc.utf8_lower(self.native))

    def to_titlecase(self) -> ArrowSeries:
        return self.with_native(pc.utf8_title(self.native))

    def zfill(self, width: int) -> ArrowSeries:
        binary_join: Incomplete = pc.binary_join_element_wise
        native = self.native
        hyphen, plus = lit("-"), lit("+")
        first_char, remaining_chars = (
            self.slice(0, 1).native,
            self.slice(1, None).native,
        )

        # Conditions
        less_than_width = pc.less(pc.utf8_length(native), lit(width))
        starts_with_hyphen = pc.equal(first_char, hyphen)
        starts_with_plus = pc.equal(first_char, plus)

        conditions = pc.make_struct(
            pc.and_(starts_with_hyphen, less_than_width),
            pc.and_(starts_with_plus, less_than_width),
            less_than_width,
        )

        # Cases
        padded_remaining_chars = pc.utf8_lpad(remaining_chars, width - 1, padding="0")

        result = pc.case_when(
            conditions,
            binary_join(
                pa.repeat(hyphen, len(native)), padded_remaining_chars, ""
            ),  # starts with hyphen and less than width
            binary_join(
                pa.repeat(plus, len(native)), padded_remaining_chars, ""
            ),  # starts with plus and less than width
            pc.utf8_lpad(native, width=width, padding="0"),  # less than width
            native,
        )
        return self.with_native(result)

    def pad_start(self, length: int, fill_char: str) -> ArrowSeries:
        return self.with_native(
            pc.utf8_lpad(self.native, width=length, padding=fill_char)
        )

    def pad_end(self, length: int, fill_char: str) -> ArrowSeries:
        return self.with_native(
            pc.utf8_rpad(self.native, width=length, padding=fill_char)
        )
