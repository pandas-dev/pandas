from __future__ import annotations

from functools import partial
import re
from typing import (
    TYPE_CHECKING,
    Callable,
    Union,
)
import warnings

import numpy as np

from pandas._libs import (
    lib,
    missing as libmissing,
)
from pandas.compat import pa_version_under7p0
from pandas.util._exceptions import find_stack_level

from pandas.core.dtypes.common import (
    is_bool_dtype,
    is_integer_dtype,
    is_object_dtype,
    is_scalar,
    is_string_dtype,
    pandas_dtype,
)
from pandas.core.dtypes.missing import isna

from pandas.core.arrays._arrow_string_mixins import ArrowStringArrayMixin
from pandas.core.arrays.arrow import ArrowExtensionArray
from pandas.core.arrays.boolean import BooleanDtype
from pandas.core.arrays.integer import Int64Dtype
from pandas.core.arrays.numeric import NumericDtype
from pandas.core.arrays.string_ import (
    BaseStringArray,
    StringDtype,
)
from pandas.core.strings.object_array import ObjectStringArrayMixin

if not pa_version_under7p0:
    import pyarrow as pa
    import pyarrow.compute as pc

    from pandas.core.arrays.arrow._arrow_utils import fallback_performancewarning


if TYPE_CHECKING:
    from pandas._typing import (
        Dtype,
        Scalar,
        npt,
    )


ArrowStringScalarOrNAT = Union[str, libmissing.NAType]


def _chk_pyarrow_available() -> None:
    if pa_version_under7p0:
        msg = "pyarrow>=7.0.0 is required for PyArrow backed ArrowExtensionArray."
        raise ImportError(msg)


# TODO: Inherit directly from BaseStringArrayMethods. Currently we inherit from
# ObjectStringArrayMixin because we want to have the object-dtype based methods as
# fallback for the ones that pyarrow doesn't yet support


class ArrowStringArray(ObjectStringArrayMixin, ArrowExtensionArray, BaseStringArray):
    """
    Extension array for string data in a ``pyarrow.ChunkedArray``.

    .. versionadded:: 1.2.0

    .. warning::

       ArrowStringArray is considered experimental. The implementation and
       parts of the API may change without warning.

    Parameters
    ----------
    values : pyarrow.Array or pyarrow.ChunkedArray
        The array of data.

    Attributes
    ----------
    None

    Methods
    -------
    None

    See Also
    --------
    :func:`pandas.array`
        The recommended function for creating a ArrowStringArray.
    Series.str
        The string methods are available on Series backed by
        a ArrowStringArray.

    Notes
    -----
    ArrowStringArray returns a BooleanArray for comparison methods.

    Examples
    --------
    >>> pd.array(['This is', 'some text', None, 'data.'], dtype="string[pyarrow]")
    <ArrowStringArray>
    ['This is', 'some text', <NA>, 'data.']
    Length: 4, dtype: string
    """

    # error: Incompatible types in assignment (expression has type "StringDtype",
    # base class "ArrowExtensionArray" defined the type as "ArrowDtype")
    _dtype: StringDtype  # type: ignore[assignment]
    _result_converter = lambda result: BooleanDtype().__from_arrow__(result)
    _storage = "pyarrow"

    def __init__(self, values) -> None:
        super().__init__(values)
        self._dtype = StringDtype(storage=self._storage)

        if not pa.types.is_string(self._pa_array.type) and not (
            pa.types.is_dictionary(self._pa_array.type)
            and pa.types.is_string(self._pa_array.type.value_type)
        ):
            raise ValueError(
                "ArrowStringArray requires a PyArrow (chunked) array of string type"
            )

    def __len__(self) -> int:
        """
        Length of this array.

        Returns
        -------
        length : int
        """
        return len(self._pa_array)

    @classmethod
    def _from_sequence(cls, scalars, dtype: Dtype | None = None, copy: bool = False):
        from pandas.core.arrays.masked import BaseMaskedArray

        _chk_pyarrow_available()

        if dtype and not (isinstance(dtype, str) and dtype == "string"):
            dtype = pandas_dtype(dtype)
            assert isinstance(dtype, StringDtype) and dtype.storage in (
                "pyarrow",
                "pyarrow_numpy",
            )

        if isinstance(scalars, BaseMaskedArray):
            # avoid costly conversion to object dtype in ensure_string_array and
            # numerical issues with Float32Dtype
            na_values = scalars._mask
            result = scalars._data
            result = lib.ensure_string_array(result, copy=copy, convert_na_value=False)
            return cls(pa.array(result, mask=na_values, type=pa.string()))
        elif isinstance(scalars, (pa.Array, pa.ChunkedArray)):
            return cls(pc.cast(scalars, pa.string()))

        # convert non-na-likes to str
        result = lib.ensure_string_array(scalars, copy=copy)
        return cls(pa.array(result, type=pa.string(), from_pandas=True))

    @classmethod
    def _from_sequence_of_strings(
        cls, strings, dtype: Dtype | None = None, copy: bool = False
    ):
        return cls._from_sequence(strings, dtype=dtype, copy=copy)

    @property
    def dtype(self) -> StringDtype:  # type: ignore[override]
        """
        An instance of 'string[pyarrow]'.
        """
        return self._dtype

    def insert(self, loc: int, item) -> ArrowStringArray:
        if not isinstance(item, str) and item is not libmissing.NA:
            raise TypeError("Scalar must be NA or str")
        return super().insert(loc, item)

    def _maybe_convert_setitem_value(self, value):
        """Maybe convert value to be pyarrow compatible."""
        if is_scalar(value):
            if isna(value):
                value = None
            elif not isinstance(value, str):
                raise TypeError("Scalar must be NA or str")
        else:
            value = np.array(value, dtype=object, copy=True)
            value[isna(value)] = None
            for v in value:
                if not (v is None or isinstance(v, str)):
                    raise TypeError("Scalar must be NA or str")
        return super()._maybe_convert_setitem_value(value)

    def isin(self, values) -> npt.NDArray[np.bool_]:
        value_set = [
            pa_scalar.as_py()
            for pa_scalar in [pa.scalar(value, from_pandas=True) for value in values]
            if pa_scalar.type in (pa.string(), pa.null())
        ]

        # short-circuit to return all False array.
        if not len(value_set):
            return np.zeros(len(self), dtype=bool)

        result = pc.is_in(self._pa_array, value_set=pa.array(value_set))
        # pyarrow 2.0.0 returned nulls, so we explicily specify dtype to convert nulls
        # to False
        return np.array(result, dtype=np.bool_)

    def astype(self, dtype, copy: bool = True):
        dtype = pandas_dtype(dtype)

        if dtype == self.dtype:
            if copy:
                return self.copy()
            return self
        elif isinstance(dtype, NumericDtype):
            data = self._pa_array.cast(pa.from_numpy_dtype(dtype.numpy_dtype))
            return dtype.__from_arrow__(data)
        elif isinstance(dtype, np.dtype) and np.issubdtype(dtype, np.floating):
            return self.to_numpy(dtype=dtype, na_value=np.nan)

        return super().astype(dtype, copy=copy)

    @property
    def _data(self):
        # dask accesses ._data directlys
        warnings.warn(
            f"{type(self).__name__}._data is a deprecated and will be removed "
            "in a future version, use ._pa_array instead",
            FutureWarning,
            stacklevel=find_stack_level(),
        )
        return self._pa_array

    # ------------------------------------------------------------------------
    # String methods interface

    # error: Incompatible types in assignment (expression has type "NAType",
    # base class "ObjectStringArrayMixin" defined the type as "float")
    _str_na_value = libmissing.NA  # type: ignore[assignment]

    def _str_map(
        self, f, na_value=None, dtype: Dtype | None = None, convert: bool = True
    ):
        # TODO: de-duplicate with StringArray method. This method is moreless copy and
        # paste.

        from pandas.arrays import (
            BooleanArray,
            IntegerArray,
        )

        if dtype is None:
            dtype = self.dtype
        if na_value is None:
            na_value = self.dtype.na_value

        mask = isna(self)
        arr = np.asarray(self)

        if is_integer_dtype(dtype) or is_bool_dtype(dtype):
            constructor: type[IntegerArray] | type[BooleanArray]
            if is_integer_dtype(dtype):
                constructor = IntegerArray
            else:
                constructor = BooleanArray

            na_value_is_na = isna(na_value)
            if na_value_is_na:
                na_value = 1
            result = lib.map_infer_mask(
                arr,
                f,
                mask.view("uint8"),
                convert=False,
                na_value=na_value,
                # error: Argument 1 to "dtype" has incompatible type
                # "Union[ExtensionDtype, str, dtype[Any], Type[object]]"; expected
                # "Type[object]"
                dtype=np.dtype(dtype),  # type: ignore[arg-type]
            )

            if not na_value_is_na:
                mask[:] = False

            return constructor(result, mask)

        elif is_string_dtype(dtype) and not is_object_dtype(dtype):
            # i.e. StringDtype
            result = lib.map_infer_mask(
                arr, f, mask.view("uint8"), convert=False, na_value=na_value
            )
            result = pa.array(result, mask=mask, type=pa.string(), from_pandas=True)
            return type(self)(result)
        else:
            # This is when the result type is object. We reach this when
            # -> We know the result type is truly object (e.g. .encode returns bytes
            #    or .findall returns a list).
            # -> We don't know the result type. E.g. `.get` can return anything.
            return lib.map_infer_mask(arr, f, mask.view("uint8"))

    def _str_contains(
        self, pat, case: bool = True, flags: int = 0, na=np.nan, regex: bool = True
    ):
        if flags:
            fallback_performancewarning()
            return super()._str_contains(pat, case, flags, na, regex)

        if regex:
            result = pc.match_substring_regex(self._pa_array, pat, ignore_case=not case)
        else:
            result = pc.match_substring(self._pa_array, pat, ignore_case=not case)
        result = self._result_converter(result)
        if not isna(na):
            result[isna(result)] = bool(na)
        return result

    def _str_startswith(self, pat: str, na=None):
        result = pc.starts_with(self._pa_array, pattern=pat)
        if not isna(na):
            result = result.fill_null(na)
        result = self._result_converter_(result)
        if not isna(na):
            result[isna(result)] = bool(na)
        return result

    def _str_endswith(self, pat: str, na=None):
        result = pc.ends_with(self._pa_array, pattern=pat)
        if not isna(na):
            result = result.fill_null(na)
        result = self._result_converter(result)
        if not isna(na):
            result[isna(result)] = bool(na)
        return result

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
            fallback_performancewarning()
            return super()._str_replace(pat, repl, n, case, flags, regex)

        func = pc.replace_substring_regex if regex else pc.replace_substring
        result = func(self._pa_array, pattern=pat, replacement=repl, max_replacements=n)
        return type(self)(result)

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

    def _str_isalnum(self):
        result = pc.utf8_is_alnum(self._pa_array)
        return self._result_converter(result)

    def _str_isalpha(self):
        result = pc.utf8_is_alpha(self._pa_array)
        return self._result_converter(result)

    def _str_isdecimal(self):
        result = pc.utf8_is_decimal(self._pa_array)
        return self._result_converter(result)

    def _str_isdigit(self):
        result = pc.utf8_is_digit(self._pa_array)
        return self._result_converter(result)

    def _str_islower(self):
        result = pc.utf8_is_lower(self._pa_array)
        return self._result_converter(result)

    def _str_isnumeric(self):
        result = pc.utf8_is_numeric(self._pa_array)
        return self._result_converter(result)

    def _str_isspace(self):
        result = pc.utf8_is_space(self._pa_array)
        return self._result_converter(result)

    def _str_istitle(self):
        result = pc.utf8_is_title(self._pa_array)
        return self._result_converter(result)

    def _str_isupper(self):
        result = pc.utf8_is_upper(self._pa_array)
        return self._result_converter(result)

    def _str_len(self):
        result = pc.utf8_length(self._pa_array)
        return Int64Dtype().__from_arrow__(result)

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


class ArrowStringArrayNumpySemantics(ArrowStringArray):
    _result_converter = lambda result: result.to_numpy()
    _storage = "pyarrow_numpy"

    def __getattribute__(self, item):
        if item in ArrowStringArrayMixin.__dict__:
            return partial(getattr(ArrowStringArrayMixin, item), self)
        return super().__getattribute__(item)

    def _str_len(self):
        result = pc.utf8_length(self._pa_array)
        return result.to_numpy()

    def _str_map(
        self, f, na_value=None, dtype: Dtype | None = None, convert: bool = True
    ):
        """
        Map a callable over valid elements of the array.

        Parameters
        ----------
        f : Callable
            A function to call on each non-NA element.
        na_value : Scalar, optional
            The value to set for NA values. Might also be used for the
            fill value if the callable `f` raises an exception.
            This defaults to ``self._str_na_value`` which is ``np.nan``
            for object-dtype and Categorical and ``pd.NA`` for StringArray.
        dtype : Dtype, optional
            The dtype of the result array.
        convert : bool, default True
            Whether to call `maybe_convert_objects` on the resulting ndarray
        """
        if dtype is None:
            dtype = np.dtype("object")
        if na_value is None:
            na_value = self._str_na_value

        if not len(self):
            return np.array([], dtype=dtype)

        arr = np.asarray(self, dtype=object)
        mask = isna(arr)
        map_convert = convert and not np.all(mask)
        try:
            result = lib.map_infer_mask(arr, f, mask.view(np.uint8), map_convert)
        except (TypeError, AttributeError) as err:
            # Reraise the exception if callable `f` got wrong number of args.
            # The user may want to be warned by this, instead of getting NaN
            p_err = (
                r"((takes)|(missing)) (?(2)from \d+ to )?\d+ "
                r"(?(3)required )positional arguments?"
            )

            if len(err.args) >= 1 and re.search(p_err, err.args[0]):
                # FIXME: this should be totally avoidable
                raise err

            def g(x):
                # This type of fallback behavior can be removed once
                # we remove object-dtype .str accessor.
                try:
                    return f(x)
                except (TypeError, AttributeError):
                    return na_value

            return self._str_map(g, na_value=na_value, dtype=dtype)
        if not isinstance(result, np.ndarray):
            return result
        if na_value is not np.nan:
            np.putmask(result, mask, na_value)
            if convert and result.dtype == object:
                result = lib.maybe_convert_objects(result)
        return result

    def _str_count(self, pat: str, flags: int = 0):
        if flags:
            return super()._str_count(pat, flags)
        return pc.count_substring_regex(self._pa_array, pat).to_numpy()

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
            return super()._str_find(sub, start, end)
        return type(self)(result)

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
        return split_func(self._pa_array, pat, max_splits=n).to_numpy()

    def _str_rsplit(self, pat: str | None = None, n: int | None = -1):
        if n in {-1, 0}:
            n = None
        return pc.split_pattern(
            self._pa_array, pat, max_splits=n, reverse=True
        ).to_numpy()
