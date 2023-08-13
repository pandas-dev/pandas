from __future__ import annotations

import re
import textwrap
from typing import (
    TYPE_CHECKING,
    Callable,
    Literal,
)
import unicodedata

import numpy as np

from pandas.compat import pa_version_under7p0

from pandas.core.dtypes.missing import isna

if TYPE_CHECKING:
    from collections.abc import Sequence

    from pandas._typing import Scalar

if not pa_version_under7p0:
    import pyarrow as pa
    import pyarrow.compute as pc


class ArrowStringArrayMixin:
    def _str_count(self, pat: str, flags: int = 0):
        if flags:
            raise NotImplementedError(f"count not implemented with {flags=}")
        return type(self)(pc.count_substring_regex(self._pa_array, pat))

    def _str_pad(
        self,
        width: int,
        side: Literal["left", "right", "both"] = "left",
        fillchar: str = " ",
    ):
        if side == "left":
            pa_pad = pc.utf8_lpad
        elif side == "right":
            pa_pad = pc.utf8_rpad
        elif side == "both":
            pa_pad = pc.utf8_center
        else:
            raise ValueError(
                f"Invalid side: {side}. Side must be one of 'left', 'right', 'both'"
            )
        return type(self)(pa_pad(self._pa_array, width=width, padding=fillchar))

    def _str_contains(
        self, pat, case: bool = True, flags: int = 0, na=None, regex: bool = True
    ):
        if flags:
            raise NotImplementedError(f"contains not implemented with {flags=}")

        if regex:
            pa_contains = pc.match_substring_regex
        else:
            pa_contains = pc.match_substring
        result = pa_contains(self._pa_array, pat, ignore_case=not case)
        if not isna(na):
            result = result.fill_null(na)
        return type(self)(result)

    def _str_startswith(self, pat: str, na=None):
        result = pc.starts_with(self._pa_array, pattern=pat)
        if not isna(na):
            result = result.fill_null(na)
        return type(self)(result)

    def _str_endswith(self, pat: str, na=None):
        result = pc.ends_with(self._pa_array, pattern=pat)
        if not isna(na):
            result = result.fill_null(na)
        return type(self)(result)

    def _str_replace(
        self,
        pat: str | re.Pattern,
        repl: str | Callable,
        n: int = -1,
        case: bool = True,
        flags: int = 0,
        regex: bool = True,
    ):
        if isinstance(pat, re.Pattern) or callable(repl) or not case or flags:
            raise NotImplementedError(
                "replace is not supported with a re.Pattern, callable repl, "
                "case=False, or flags!=0"
            )

        func = pc.replace_substring_regex if regex else pc.replace_substring
        result = func(self._pa_array, pattern=pat, replacement=repl, max_replacements=n)
        return type(self)(result)

    def _str_repeat(self, repeats: int | Sequence[int]):
        if not isinstance(repeats, int):
            raise NotImplementedError(
                f"repeat is not implemented when repeats is {type(repeats).__name__}"
            )
        else:
            return type(self)(pc.binary_repeat(self._pa_array, repeats))

    def _str_match(
        self, pat: str, case: bool = True, flags: int = 0, na: Scalar | None = None
    ):
        if not pat.startswith("^"):
            pat = f"^{pat}"
        return self._str_contains(pat, case, flags, na, regex=True)

    def _str_fullmatch(
        self, pat, case: bool = True, flags: int = 0, na: Scalar | None = None
    ):
        if not pat.endswith("$") or pat.endswith("//$"):
            pat = f"{pat}$"
        return self._str_match(pat, case, flags, na)

    def _str_find(self, sub: str, start: int = 0, end: int | None = None):
        if start != 0 and end is not None:
            slices = pc.utf8_slice_codeunits(self._pa_array, start, stop=end)
            result = pc.find_substring(slices, sub)
            not_found = pc.equal(result, -1)
            offset_result = pc.add(result, end - start)
            result = pc.if_else(not_found, result, offset_result)
        elif start == 0 and end is None:
            slices = self._pa_array
            result = pc.find_substring(slices, sub)
        else:
            raise NotImplementedError(
                f"find not implemented with {sub=}, {start=}, {end=}"
            )
        return type(self)(result)

    def _str_get(self, i: int):
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
        return type(self)(result)

    def _str_join(self, sep: str):
        if pa.types.is_string(self._pa_array.type):
            result = self._apply_elementwise(list)
            result = pa.chunked_array(result, type=pa.list_(pa.string()))
        else:
            result = self._pa_array
        return type(self)(pc.binary_join(result, sep))

    def _str_partition(self, sep: str, expand: bool):
        predicate = lambda val: val.partition(sep)
        result = self._apply_elementwise(predicate)
        return type(self)(pa.chunked_array(result))

    def _str_rpartition(self, sep: str, expand: bool):
        predicate = lambda val: val.rpartition(sep)
        result = self._apply_elementwise(predicate)
        return type(self)(pa.chunked_array(result))

    def _str_slice(
        self, start: int | None = None, stop: int | None = None, step: int | None = None
    ):
        if start is None:
            start = 0
        if step is None:
            step = 1
        return type(self)(
            pc.utf8_slice_codeunits(self._pa_array, start=start, stop=stop, step=step)
        )

    def _str_slice_replace(
        self, start: int | None = None, stop: int | None = None, repl: str | None = None
    ):
        if repl is None:
            repl = ""
        if start is None:
            start = 0
        return type(self)(pc.utf8_replace_slice(self._pa_array, start, stop, repl))

    def _str_isalnum(self):
        return type(self)(pc.utf8_is_alnum(self._pa_array))

    def _str_isalpha(self):
        return type(self)(pc.utf8_is_alpha(self._pa_array))

    def _str_isdecimal(self):
        return type(self)(pc.utf8_is_decimal(self._pa_array))

    def _str_isdigit(self):
        return type(self)(pc.utf8_is_digit(self._pa_array))

    def _str_islower(self):
        return type(self)(pc.utf8_is_lower(self._pa_array))

    def _str_isnumeric(self):
        return type(self)(pc.utf8_is_numeric(self._pa_array))

    def _str_isspace(self):
        return type(self)(pc.utf8_is_space(self._pa_array))

    def _str_istitle(self):
        return type(self)(pc.utf8_is_title(self._pa_array))

    def _str_capitalize(self):
        return type(self)(pc.utf8_capitalize(self._pa_array))

    def _str_title(self):
        return type(self)(pc.utf8_title(self._pa_array))

    def _str_isupper(self):
        return type(self)(pc.utf8_is_upper(self._pa_array))

    def _str_swapcase(self):
        return type(self)(pc.utf8_swapcase(self._pa_array))

    def _str_len(self):
        return type(self)(pc.utf8_length(self._pa_array))

    def _str_lower(self):
        return type(self)(pc.utf8_lower(self._pa_array))

    def _str_upper(self):
        return type(self)(pc.utf8_upper(self._pa_array))

    def _str_strip(self, to_strip=None):
        if to_strip is None:
            result = pc.utf8_trim_whitespace(self._pa_array)
        else:
            result = pc.utf8_trim(self._pa_array, characters=to_strip)
        return type(self)(result)

    def _str_lstrip(self, to_strip=None):
        if to_strip is None:
            result = pc.utf8_ltrim_whitespace(self._pa_array)
        else:
            result = pc.utf8_ltrim(self._pa_array, characters=to_strip)
        return type(self)(result)

    def _str_rstrip(self, to_strip=None):
        if to_strip is None:
            result = pc.utf8_rtrim_whitespace(self._pa_array)
        else:
            result = pc.utf8_rtrim(self._pa_array, characters=to_strip)
        return type(self)(result)

    def _str_removeprefix(self, prefix: str):
        # TODO: Should work once https://github.com/apache/arrow/issues/14991 is fixed
        # starts_with = pc.starts_with(self._pa_array, pattern=prefix)
        # removed = pc.utf8_slice_codeunits(self._pa_array, len(prefix))
        # result = pc.if_else(starts_with, removed, self._pa_array)
        # return type(self)(result)
        predicate = lambda val: val.removeprefix(prefix)
        result = self._apply_elementwise(predicate)
        return type(self)(pa.chunked_array(result))

    def _str_removesuffix(self, suffix: str):
        ends_with = pc.ends_with(self._pa_array, pattern=suffix)
        removed = pc.utf8_slice_codeunits(self._pa_array, 0, stop=-len(suffix))
        result = pc.if_else(ends_with, removed, self._pa_array)
        return type(self)(result)

    def _str_casefold(self):
        predicate = lambda val: val.casefold()
        result = self._apply_elementwise(predicate)
        return type(self)(pa.chunked_array(result))

    def _str_encode(self, encoding: str, errors: str = "strict"):
        predicate = lambda val: val.encode(encoding, errors)
        result = self._apply_elementwise(predicate)
        return type(self)(pa.chunked_array(result))

    def _str_extract(self, pat: str, flags: int = 0, expand: bool = True):
        raise NotImplementedError(
            "str.extract not supported with pd.ArrowDtype(pa.string())."
        )

    def _str_findall(self, pat: str, flags: int = 0):
        regex = re.compile(pat, flags=flags)
        predicate = lambda val: regex.findall(val)
        result = self._apply_elementwise(predicate)
        return type(self)(pa.chunked_array(result))

    def _str_get_dummies(self, sep: str = "|"):
        split = pc.split_pattern(self._pa_array, sep)
        flattened_values = pc.list_flatten(split)
        uniques = flattened_values.unique()
        uniques_sorted = uniques.take(pa.compute.array_sort_indices(uniques))
        lengths = pc.list_value_length(split).fill_null(0).to_numpy()
        n_rows = len(self)
        n_cols = len(uniques)
        indices = pc.index_in(flattened_values, uniques_sorted).to_numpy()
        indices = indices + np.arange(n_rows).repeat(lengths) * n_cols
        dummies = np.zeros(n_rows * n_cols, dtype=np.bool_)
        dummies[indices] = True
        dummies = dummies.reshape((n_rows, n_cols))
        result = type(self)(pa.array(list(dummies)))
        return result, uniques_sorted.to_pylist()

    def _str_index(self, sub: str, start: int = 0, end: int | None = None):
        predicate = lambda val: val.index(sub, start, end)
        result = self._apply_elementwise(predicate)
        return type(self)(pa.chunked_array(result))

    def _str_rindex(self, sub: str, start: int = 0, end: int | None = None):
        predicate = lambda val: val.rindex(sub, start, end)
        result = self._apply_elementwise(predicate)
        return type(self)(pa.chunked_array(result))

    def _str_normalize(self, form: str):
        predicate = lambda val: unicodedata.normalize(form, val)
        result = self._apply_elementwise(predicate)
        return type(self)(pa.chunked_array(result))

    def _str_rfind(self, sub: str, start: int = 0, end=None):
        predicate = lambda val: val.rfind(sub, start, end)
        result = self._apply_elementwise(predicate)
        return type(self)(pa.chunked_array(result))

    def _str_split(
        self,
        pat: str | None = None,
        n: int | None = -1,
        expand: bool = False,
        regex: bool | None = None,
    ):
        if n in {-1, 0}:
            n = None
        if regex:
            split_func = pc.split_pattern_regex
        else:
            split_func = pc.split_pattern
        return type(self)(split_func(self._pa_array, pat, max_splits=n))

    def _str_rsplit(self, pat: str | None = None, n: int | None = -1):
        if n in {-1, 0}:
            n = None
        return type(self)(
            pc.split_pattern(self._pa_array, pat, max_splits=n, reverse=True)
        )

    def _str_translate(self, table: dict[int, str]):
        predicate = lambda val: val.translate(table)
        result = self._apply_elementwise(predicate)
        return type(self)(pa.chunked_array(result))

    def _str_wrap(self, width: int, **kwargs):
        kwargs["width"] = width
        tw = textwrap.TextWrapper(**kwargs)
        predicate = lambda val: "\n".join(tw.wrap(val))
        result = self._apply_elementwise(predicate)
        return type(self)(pa.chunked_array(result))
