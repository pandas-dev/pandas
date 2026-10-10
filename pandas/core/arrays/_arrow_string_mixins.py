from __future__ import annotations

from functools import partial
import re
from typing import (
    TYPE_CHECKING,
    Any,
    Literal,
    Self,
)
import unicodedata

import numpy as np

from pandas._libs import lib
from pandas.compat import (
    HAS_PYARROW,
    pa_version_under17p0,
    pa_version_under21p0,
)

if HAS_PYARROW:
    import pyarrow as pa
    import pyarrow.compute as pc

from pandas.core.arrays._re2 import (
    RE2Pattern,
    translate,
    translate_template,
)

if TYPE_CHECKING:
    from collections.abc import Callable

    from pandas._typing import Scalar


class ArrowStringArrayMixin:
    _pa_array: pa.ChunkedArray

    def __init__(self, *args, **kwargs) -> None:
        raise NotImplementedError

    def _from_pyarrow_array(self, pa_array: pa.Array | pa.ChunkedArray) -> Self:
        raise NotImplementedError

    def _convert_bool_result(self, result, na=lib.no_default, method_name=None):
        # Convert a bool-dtype result to the appropriate result type
        raise NotImplementedError

    def _convert_int_result(self, result):
        # Convert an integer-dtype result to the appropriate result type
        raise NotImplementedError

    def _apply_elementwise(self, func: Callable[..., Any]) -> list[list[Any]]:
        raise NotImplementedError

    def _str_re_fallback(self, method: str, *args, **kwargs):
        # Evaluate the ``_str_<method>`` regex method with Python's ``re`` for a
        #  pattern or data on which pyarrow's RE2 kernels would differ
        raise NotImplementedError

    def _to_re2(
        self,
        pat: str | re.Pattern[str],
        case: bool = True,
        flags: int = 0,
        named_groups: bool = False,
    ) -> RE2Pattern | None:
        """
        Translate `pat` for pyarrow's regex kernels.

        Returns None if pyarrow cannot evaluate `pat` as ``re`` does on the
        values of this array; the caller then falls back to ``re``.

        Raises
        ------
        re.error
            If `pat` is not a valid Python pattern.
        """
        if isinstance(pat, re.Pattern):
            if flags or not case:
                # whether these combine with the flags of `pat` is decided (or
                #  rejected) by the fallback
                return None
            pattern, flags = pat.pattern, pat.flags
        else:
            pattern = pat
            if not case:
                flags |= re.IGNORECASE
        re2 = translate(pattern, flags, named_groups)
        if re2 is None or not re2.matches_python_on(self._pa_array):
            return None
        return re2

    def _str_len(self):
        result = pc.utf8_length(self._pa_array)
        return self._convert_int_result(result)

    def _str_lower(self) -> Self:
        return self._from_pyarrow_array(pc.utf8_lower(self._pa_array))

    def _str_upper(self) -> Self:
        return self._from_pyarrow_array(pc.utf8_upper(self._pa_array))

    def _str_strip(self, to_strip=None) -> Self:
        if to_strip is None:
            result = pc.utf8_trim_whitespace(self._pa_array)
        else:
            result = pc.utf8_trim(self._pa_array, characters=to_strip)
        return self._from_pyarrow_array(result)

    def _str_lstrip(self, to_strip=None) -> Self:
        if to_strip is None:
            result = pc.utf8_ltrim_whitespace(self._pa_array)
        else:
            result = pc.utf8_ltrim(self._pa_array, characters=to_strip)
        return self._from_pyarrow_array(result)

    def _str_rstrip(self, to_strip=None) -> Self:
        if to_strip is None:
            result = pc.utf8_rtrim_whitespace(self._pa_array)
        else:
            result = pc.utf8_rtrim(self._pa_array, characters=to_strip)
        return self._from_pyarrow_array(result)

    def _str_pad(
        self,
        width: int,
        side: Literal["left", "right", "both"] = "left",
        fillchar: str = " ",
    ) -> Self:
        if width < 0:
            width = 0
        if side == "left":
            pa_pad = pc.utf8_lpad
        elif side == "right":
            pa_pad = pc.utf8_rpad
        elif side == "both":
            if pa_version_under17p0:
                # GH#59624 fall back to object dtype
                from pandas import array

                obj_arr = self.astype(object, copy=False)  # type: ignore[attr-defined]
                obj = array(obj_arr, dtype=object)
                result = obj._str_pad(width, side, fillchar)  # type: ignore[attr-defined]
                return type(self)._from_sequence(result, dtype=self.dtype)  # type: ignore[attr-defined]
            else:
                # GH#54792
                # https://github.com/apache/arrow/issues/15053#issuecomment-2317032347
                lean_left = (width % 2) == 0
                pa_pad = partial(pc.utf8_center, lean_left_on_odd_padding=lean_left)
        else:
            raise ValueError(
                f"Invalid side: {side}. Side must be one of 'left', 'right', 'both'"
            )
        return self._from_pyarrow_array(
            pa_pad(self._pa_array, width=width, padding=fillchar)
        )

    def _str_zfill(self, width: int) -> Self:
        if pa_version_under21p0:
            predicate = lambda val: val.zfill(width)
            result = self._apply_elementwise(predicate)
            return self._from_pyarrow_array(
                pa.chunked_array(result, type=self._pa_array.type)
            )
        # pc.utf8_zfill raises on a negative width, while str.zfill returns the
        # string unchanged -> clamping to zero to get that behaviour (GH#69486)
        width = max(width, 0)
        return self._from_pyarrow_array(pc.utf8_zfill(self._pa_array, width))

    def _str_normalize(self, form: Literal["NFC", "NFD", "NFKC", "NFKD"]) -> Self:
        if form not in ("NFC", "NFD", "NFKC", "NFKD"):
            raise ValueError("invalid normalization form")
        if form in ("NFC", "NFKC"):
            # GH#64359 pc.utf8_normalize only decomposes; it skips the canonical
            #  composition step, so for the composing forms it returns decomposed
            #  output. Fall back to unicodedata for these.
            predicate = lambda val: unicodedata.normalize(form, val)
            result = self._apply_elementwise(predicate)
            return self._from_pyarrow_array(
                pa.chunked_array(result, type=self._pa_array.type)
            )
        return self._from_pyarrow_array(pc.utf8_normalize(self._pa_array, form=form))

    def _str_get(self, i: int) -> Self:
        lengths = pc.utf8_length(self._pa_array)
        if i >= 0:
            out_of_bounds = pc.greater_equal(i, lengths)
            start = i
            stop = i + 1
            step = 1
        else:
            out_of_bounds = pc.greater(-i, lengths)
            start = i
            stop = i - 1
            step = -1
        not_out_of_bounds = pc.invert(out_of_bounds.fill_null(True))
        selected = pc.utf8_slice_codeunits(
            self._pa_array, start=start, stop=stop, step=step
        )
        null_value = pa.scalar(None, type=self._pa_array.type)
        result = pc.if_else(not_out_of_bounds, selected, null_value)
        return self._from_pyarrow_array(result)

    def _str_slice(
        self, start: int | None = None, stop: int | None = None, step: int | None = None
    ) -> Self:
        if start is None:
            if step is not None and step < 0:
                # GH#59710
                start = -1
            else:
                start = 0
        if step is None:
            step = 1
        return self._from_pyarrow_array(
            pc.utf8_slice_codeunits(self._pa_array, start=start, stop=stop, step=step)
        )

    def _str_getitem(self, key: slice | int) -> Self:
        if isinstance(key, slice):
            return self._str_slice(start=key.start, stop=key.stop, step=key.step)
        else:
            return self._str_get(key)

    def _str_slice_replace(
        self, start: int | None = None, stop: int | None = None, repl: str | None = None
    ) -> Self:
        if repl is None:
            repl = ""
        if start is None:
            start = 0
        if stop is None:
            stop = np.iinfo(np.int64).max
        return self._from_pyarrow_array(
            pc.utf8_replace_slice(self._pa_array, start, stop, repl)
        )

    def _str_replace(
        self,
        pat: str | re.Pattern[str],
        repl: str | Callable[..., Any],
        n: int = -1,
        case: bool = True,
        flags: int = 0,
        regex: bool = True,
    ) -> Self:
        # https://github.com/apache/arrow/issues/39149
        # GH 56404, unexpected behavior with negative max_replacements with pyarrow.
        pa_max_replacements = None if n < 0 else n

        if regex and isinstance(repl, str):
            re2 = self._to_re2(pat, case, flags)
            # pyarrow steps over empty matches differently from re, and with
            #  max_replacements it re-evaluates assertions without the text
            #  around the match
            if (
                re2 is not None
                and not re2.nullable
                and not re2.nullable_loop
                and not (n >= 0 and re2.has_assertions)
            ):
                compiled = re.compile(pat, flags if case else flags | re.IGNORECASE)
                rewrite = translate_template(compiled, repl)
                if rewrite is not None:
                    result = pc.replace_substring_regex(
                        self._pa_array,
                        pattern=re2.pattern,
                        replacement=rewrite[0],
                        max_replacements=pa_max_replacements,
                    )
                    return self._from_pyarrow_array(result)
        elif (
            not regex
            and case
            and not flags
            and isinstance(pat, str)
            and isinstance(repl, str)
        ):
            if pat == "":
                # pyarrow hangs for empty patterns
                # (https://github.com/apache/arrow/issues/39149)
                func = lambda val: val.replace(pat, repl, n)
                result = self._apply_elementwise(func)
                return self._from_pyarrow_array(
                    pa.chunked_array(result, type=self._pa_array.type)
                )
            result = pc.replace_substring(
                self._pa_array,
                pattern=pat,
                replacement=repl,
                max_replacements=pa_max_replacements,
            )
            return self._from_pyarrow_array(result)
        return self._str_re_fallback("replace", pat, repl, n, case, flags, regex)

    def _str_capitalize(self) -> Self:
        return self._from_pyarrow_array(pc.utf8_capitalize(self._pa_array))

    def _str_title(self) -> Self:
        return self._from_pyarrow_array(pc.utf8_title(self._pa_array))

    def _str_swapcase(self) -> Self:
        return self._from_pyarrow_array(pc.utf8_swapcase(self._pa_array))

    def _str_removeprefix(self, prefix: str):
        if prefix == "":
            return self._from_pyarrow_array(self._pa_array)
        starts_with = pc.starts_with(self._pa_array, pattern=prefix)
        removed = pc.utf8_slice_codeunits(self._pa_array, len(prefix))
        result = pc.if_else(starts_with, removed, self._pa_array)
        return self._from_pyarrow_array(result)

    def _str_removesuffix(self, suffix: str):
        if suffix == "":
            return self._from_pyarrow_array(self._pa_array)
        ends_with = pc.ends_with(self._pa_array, pattern=suffix)
        removed = pc.utf8_slice_codeunits(self._pa_array, 0, stop=-len(suffix))
        result = pc.if_else(ends_with, removed, self._pa_array)
        return self._from_pyarrow_array(result)

    def _str_startswith(
        self, pat: str | tuple[str, ...], na: Scalar | lib.NoDefault = lib.no_default
    ):
        if isinstance(pat, str):
            result = pc.starts_with(self._pa_array, pattern=pat)
        elif len(pat) == 0:
            # For empty tuple we return null for missing values and False
            #  for valid values.
            result = pc.if_else(pc.is_null(self._pa_array), None, False)
        else:
            result = pc.starts_with(self._pa_array, pattern=pat[0])

            for p in pat[1:]:
                result = pc.or_(result, pc.starts_with(self._pa_array, pattern=p))
        return self._convert_bool_result(result, na=na, method_name="startswith")

    def _str_endswith(
        self, pat: str | tuple[str, ...], na: Scalar | lib.NoDefault = lib.no_default
    ):
        if isinstance(pat, str):
            result = pc.ends_with(self._pa_array, pattern=pat)
        elif len(pat) == 0:
            # For empty tuple we return null for missing values and False
            #  for valid values.
            result = pc.if_else(pc.is_null(self._pa_array), None, False)
        else:
            result = pc.ends_with(self._pa_array, pattern=pat[0])

            for p in pat[1:]:
                result = pc.or_(result, pc.ends_with(self._pa_array, pattern=p))
        return self._convert_bool_result(result, na=na, method_name="endswith")

    def _str_isalnum(self):
        result = pc.utf8_is_alnum(self._pa_array)
        return self._convert_bool_result(result)

    def _str_isalpha(self):
        result = pc.utf8_is_alpha(self._pa_array)
        return self._convert_bool_result(result)

    def _str_isascii(self):
        result = pc.string_is_ascii(self._pa_array)
        return self._convert_bool_result(result)

    def _str_isdecimal(self):
        result = pc.utf8_is_decimal(self._pa_array)
        return self._convert_bool_result(result)

    def _str_isdigit(self):
        if pa_version_under21p0:
            # https://github.com/pandas-dev/pandas/issues/61466
            res_list = self._apply_elementwise(str.isdigit)
            return self._convert_bool_result(
                pa.chunked_array(res_list, type=pa.bool_())
            )
        result = pc.utf8_is_digit(self._pa_array)
        return self._convert_bool_result(result)

    def _str_islower(self):
        result = pc.utf8_is_lower(self._pa_array)
        return self._convert_bool_result(result)

    def _str_isnumeric(self):
        result = pc.utf8_is_numeric(self._pa_array)
        return self._convert_bool_result(result)

    def _str_isspace(self):
        result = pc.utf8_is_space(self._pa_array)
        return self._convert_bool_result(result)

    def _str_istitle(self):
        result = pc.utf8_is_title(self._pa_array)
        return self._convert_bool_result(result)

    def _str_isupper(self):
        result = pc.utf8_is_upper(self._pa_array)
        return self._convert_bool_result(result)

    def _str_contains(
        self,
        pat,
        case: bool = True,
        flags: int = 0,
        na: Scalar | lib.NoDefault = lib.no_default,
        regex: bool = True,
    ):
        if regex:
            re2 = self._to_re2(pat, case, flags)
            if re2 is None:
                return self._str_re_fallback("contains", pat, case, flags, na, regex)
            result = pc.match_substring_regex(self._pa_array, re2.pattern)
        elif flags or isinstance(pat, re.Pattern):
            return self._str_re_fallback("contains", pat, case, flags, na, regex)
        else:
            result = pc.match_substring(self._pa_array, pat, ignore_case=not case)
        return self._convert_bool_result(result, na=na, method_name="contains")

    def _str_match(
        self,
        pat: str | re.Pattern[str],
        case: bool = True,
        flags: int = 0,
        na: Scalar | lib.NoDefault = lib.no_default,
    ):
        re2 = self._to_re2(pat, case, flags)
        if re2 is None:
            return self._str_re_fallback("match", pat, case, flags, na)
        result = pc.match_substring_regex(self._pa_array, rf"\A(?:{re2.pattern})")
        return self._convert_bool_result(result, na=na, method_name="match")

    def _str_fullmatch(
        self,
        pat: str | re.Pattern[str],
        case: bool = True,
        flags: int = 0,
        na: Scalar | lib.NoDefault = lib.no_default,
    ):
        re2 = self._to_re2(pat, case, flags)
        if re2 is None:
            return self._str_re_fallback("fullmatch", pat, case, flags, na)
        result = pc.match_substring_regex(self._pa_array, rf"\A(?:{re2.pattern})\z")
        return self._convert_bool_result(result, na=na, method_name="fullmatch")

    def _str_count(self, pat: str | re.Pattern[str], flags: int = 0):
        re2 = self._to_re2(pat, flags=flags)
        # pyarrow restarts the search on the remainder of the string after each
        #  match, which loses the text before it, and steps over empty matches
        #  by bytes rather than characters
        if re2 is None or re2.nullable or re2.nullable_loop or re2.start_assertions:
            return self._str_re_fallback("count", pat, flags)
        result = pc.count_substring_regex(self._pa_array, re2.pattern)
        return self._convert_int_result(result)

    def _str_find(self, sub: str, start: int = 0, end: int | None = None):
        # min_count=0 so that an empty or all-null array reports True instead of
        #  null, keeping it on the pyarrow path below
        if not pc.all(pc.string_is_ascii(self._pa_array), min_count=0).as_py():
            # GH#64123 - pc.find_substring returns byte offsets instead of
            # character offsets for multi-byte UTF-8 characters, so we fall back
            # to Python str.find which correctly returns character offsets.
            res_list = self._apply_elementwise(lambda val: val.find(sub, start, end))
            return self._convert_int_result(pa.chunked_array(res_list, type=pa.int64()))

        if (start == 0 or start is None) and end is None:
            result = pc.find_substring(self._pa_array, sub)
        else:
            if sub == "":
                # GH#56792
                res_list = self._apply_elementwise(
                    lambda val: val.find(sub, start, end)
                )
                return self._convert_int_result(
                    pa.chunked_array(res_list, type=pa.int64())
                )
            if start is None:
                start_offset = 0
                start = 0
            elif start < 0:
                start_offset = pc.add(start, pc.utf8_length(self._pa_array))
                start_offset = pc.if_else(pc.less(start_offset, 0), 0, start_offset)
            else:
                start_offset = start
            slices = pc.utf8_slice_codeunits(self._pa_array, start, stop=end)
            result = pc.find_substring(slices, sub)
            found = pc.not_equal(result, pa.scalar(-1, type=result.type))
            offset_result = pc.add(result, start_offset)
            result = pc.if_else(found, offset_result, -1)
        result = result.cast(pa.int64())
        return self._convert_int_result(result)

    def _str_partition_expand(self, sep: str) -> pa.ChunkedArray:
        """
        Split each string on the first occurrence of ``sep``.

        Returns a ``list<string>`` array holding one three-element list per row
        -- the part before the separator, the separator, and the part after --
        which ``StringMethods._wrap_result`` expands into three columns. Rows
        without ``sep`` get two empty strings, matching ``str.partition``.

        The caller wraps this in an :class:`ArrowExtensionArray`; the rows are
        lists, so it is not of the calling array's own type.
        """
        if not sep:
            # pyarrow reports this as "Empty separator"; keep str.partition's
            #  wording so every dtype raises the same way
            raise ValueError("empty separator")

        str_type = self._pa_array.type
        chunks = [
            self._partition_chunk(chunk, sep, str_type)
            for chunk in self._pa_array.chunks
        ]
        return pa.chunked_array(chunks, type=pa.list_(str_type))

    @staticmethod
    def _partition_chunk(chunk: pa.Array, sep: str, str_type: pa.DataType) -> pa.Array:
        """
        Build the ``list<string>`` rows for one chunk of :meth:`_str_partition_expand`.

        Working a chunk at a time keeps the concatenation below within the
        offset width of ``str_type``, which matters for columns near the 2 GiB
        limit of 32-bit ``string``.
        """
        # max_splits=1 gives [before] when sep is absent, else [before, after];
        #  padding to a fixed two elements makes the tail null in the first case
        split = pc.split_pattern(chunk, sep, max_splits=1)
        pieces = pc.list_slice(split, 0, 2, return_fixed_size_list=True)

        before = pc.list_element(pieces, 0)
        tail = pc.list_element(pieces, 1)
        after = pc.fill_null(tail, pa.scalar("", type=str_type))
        middle = pc.if_else(
            pc.is_valid(tail),
            pa.scalar(sep, type=str_type),
            pa.scalar("", type=str_type),
        )

        n = len(chunk)
        # Interleave the three columns into [before[0], sep[0], after[0],
        #  before[1], ...]. Arrow has no kernel that zips arrays into rows, so
        #  take from their concatenation. Both ranges below are plain index
        #  arithmetic: NumPy builds them in one strided pass, while pyarrow has
        #  no arange at all before 21.0 and would need four kernels for the
        #  interleave, so they stay NumPy at every supported version.
        values = pa.concat_arrays([before, middle, after])
        indices = np.arange(3 * n, dtype=np.int64).reshape(3, n).T.reshape(-1)
        values = values.take(pa.array(indices))
        # int32 rather than int64 offsets, but built wide so that a chunk with
        #  more than 2**31 / 3 rows raises instead of wrapping around
        offsets = pa.array(np.arange(0, 3 * n + 1, 3, dtype=np.int64), type=pa.int32())
        # a null string partitions to a null row, not to a row of nulls
        mask = pc.is_null(before) if before.null_count else None
        return pa.ListArray.from_arrays(offsets, values, mask=mask)
