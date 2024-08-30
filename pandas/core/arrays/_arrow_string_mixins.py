from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Literal,
)

import numpy as np

from pandas.compat import pa_version_under10p1

from pandas.core.dtypes.missing import isna

if not pa_version_under10p1:
    import pyarrow as pa
    import pyarrow.compute as pc

if TYPE_CHECKING:
    from collections.abc import Sized

    from pandas._typing import (
        Scalar,
        Self,
    )


class ArrowStringArrayMixin:
    _pa_array: Sized

    def __init__(self, *args, **kwargs) -> None:
        raise NotImplementedError

    def _convert_bool_result(self, result):
        # Convert a bool-dtype result to the appropriate result type
        raise NotImplementedError

    def _convert_int_result(self, result):
        # Convert an integer-dtype result to the appropriate result type
        raise NotImplementedError

    def _str_pad(
        self,
        width: int,
        side: Literal["left", "right", "both"] = "left",
        fillchar: str = " ",
    ) -> Self:
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
        null_value = pa.scalar(
            None,
            type=self._pa_array.type,  # type: ignore[attr-defined]
        )
        result = pc.if_else(not_out_of_bounds, selected, null_value)
        return type(self)(result)

    def _str_slice_replace(
        self, start: int | None = None, stop: int | None = None, repl: str | None = None
    ) -> Self:
        if repl is None:
            repl = ""
        if start is None:
            start = 0
        if stop is None:
            stop = np.iinfo(np.int64).max
        return type(self)(pc.utf8_replace_slice(self._pa_array, start, stop, repl))

    def _str_capitalize(self) -> Self:
        return type(self)(pc.utf8_capitalize(self._pa_array))

    def _str_title(self) -> Self:
        return type(self)(pc.utf8_title(self._pa_array))

    def _str_swapcase(self) -> Self:
        return type(self)(pc.utf8_swapcase(self._pa_array))

    def _str_removesuffix(self, suffix: str):
        ends_with = pc.ends_with(self._pa_array, pattern=suffix)
        removed = pc.utf8_slice_codeunits(self._pa_array, 0, stop=-len(suffix))
        result = pc.if_else(ends_with, removed, self._pa_array)
        return type(self)(result)

    def _str_startswith(self, pat: str | tuple[str, ...], na: Scalar | None = None):
        if isinstance(pat, str):
            result = pc.starts_with(self._pa_array, pattern=pat)
        else:
            if len(pat) == 0:
                # For empty tuple we return null for missing values and False
                #  for valid values.
                result = pc.if_else(pc.is_null(self._pa_array), None, False)
            else:
                result = pc.starts_with(self._pa_array, pattern=pat[0])

                for p in pat[1:]:
                    result = pc.or_(result, pc.starts_with(self._pa_array, pattern=p))
        if not isna(na):  # pyright: ignore [reportGeneralTypeIssues]
            result = result.fill_null(na)
        return self._convert_bool_result(result)

    def _str_endswith(self, pat: str | tuple[str, ...], na: Scalar | None = None):
        if isinstance(pat, str):
            result = pc.ends_with(self._pa_array, pattern=pat)
        else:
            if len(pat) == 0:
                # For empty tuple we return null for missing values and False
                #  for valid values.
                result = pc.if_else(pc.is_null(self._pa_array), None, False)
            else:
                result = pc.ends_with(self._pa_array, pattern=pat[0])

                for p in pat[1:]:
                    result = pc.or_(result, pc.ends_with(self._pa_array, pattern=p))
        if not isna(na):  # pyright: ignore [reportGeneralTypeIssues]
            result = result.fill_null(na)
        return self._convert_bool_result(result)
