import re

import numpy as np

import pandas._libs.lib as lib

from pandas.core.dtypes.generic import ABCSeries
from pandas.core.dtypes.missing import isna

from pandas.core.accessor import CachedAccessor
from pandas.core.arrays.numpy_ import PandasArray
from pandas.core.strings.base import BaseStringArrayMethods


class ObjectArrayMethods(BaseStringArrayMethods):
    def _map(self, f, na_mask=True, na_value=np.nan, dtype=np.dtype(object)):
        # n.b.: na_mask is the default.
        # need to figure out when it was false, maybe split
        arr = self._array  # object-dtype ndarray.

        if not len(arr):
            return np.ndarray(0, dtype=dtype)
        if na_value is None:
            na_value = np.nan

        if isinstance(arr, ABCSeries):
            arr = arr._values  # TODO: extract_array?
        if not isinstance(arr, np.ndarray):
            arr = np.asarray(arr, dtype=object)
        if na_mask:
            mask = isna(arr)
            convert = not np.all(mask)
            try:
                result = lib.map_infer_mask(arr, f, mask.view(np.uint8), convert)
            except (TypeError, AttributeError) as e:
                # Reraise the exception if callable `f` got wrong number of args.
                # The user may want to be warned by this, instead of getting NaN
                p_err = (
                    r"((takes)|(missing)) (?(2)from \d+ to )?\d+ "
                    r"(?(3)required )positional arguments?"
                )

                if len(e.args) >= 1 and re.search(p_err, e.args[0]):
                    # FIXME: this should be totally avoidable
                    raise e

                def g(x):
                    try:
                        return f(x)
                    except (TypeError, AttributeError):
                        return na_value

                return self._map_object(g, dtype=dtype)
            if na_value is not np.nan:
                np.putmask(result, mask, na_value)
                if result.dtype == object:
                    result = lib.maybe_convert_objects(result)
            return result
        else:
            return lib.map_infer(arr, f)

    def cat(self, others=None, sep=None, na_rep=None, join="left"):
        from pandas import Index, Series, concat

        if isinstance(others, str):
            raise ValueError("Did you mean to supply a `sep` keyword?")
        if sep is None:
            sep = ""

        if isinstance(self._orig, ABCIndexClass):
            data = Series(self._orig, index=self._orig)
        else:  # Series
            data = self._orig

        # concatenate Series/Index with itself if no "others"
        if others is None:
            data = ensure_object(data)
            na_mask = isna(data)
            if na_rep is None and na_mask.any():
                data = data[~na_mask]
            elif na_rep is not None and na_mask.any():
                data = np.where(na_mask, na_rep, data)
            return sep.join(data)

        try:
            # turn anything in "others" into lists of Series
            others = self._get_series_list(others)
        except ValueError as err:  # do not catch TypeError raised by _get_series_list
            raise ValueError(
                "If `others` contains arrays or lists (or other "
                "list-likes without an index), these must all be "
                "of the same length as the calling Series/Index."
            ) from err

        # align if required
        if any(not data.index.equals(x.index) for x in others):
            # Need to add keys for uniqueness in case of duplicate columns
            others = concat(
                others,
                axis=1,
                join=(join if join == "inner" else "outer"),
                keys=range(len(others)),
                sort=False,
                copy=False,
            )
            data, others = data.align(others, join=join)
            others = [others[x] for x in others]  # again list of Series

        all_cols = [ensure_object(x) for x in [data] + others]
        na_masks = np.array([isna(x) for x in all_cols])
        union_mask = np.logical_or.reduce(na_masks, axis=0)

        if na_rep is None and union_mask.any():
            # no na_rep means NaNs for all rows where any column has a NaN
            # only necessary if there are actually any NaNs
            result = np.empty(len(data), dtype=object)
            np.putmask(result, union_mask, np.nan)

            not_masked = ~union_mask
            result[not_masked] = cat_safe([x[not_masked] for x in all_cols], sep)
        elif na_rep is not None and union_mask.any():
            # fill NaNs with na_rep in case there are actually any NaNs
            all_cols = [
                np.where(nm, na_rep, col) for nm, col in zip(na_masks, all_cols)
            ]
            result = cat_safe(all_cols, sep)
        else:
            # no NaNs - can just concatenate
            result = cat_safe(all_cols, sep)

        if isinstance(self._orig, ABCIndexClass):
            # add dtype for case that result is all-NA
            result = Index(result, dtype=object, name=self._orig.name)
        else:  # Series
            if is_categorical_dtype(self._orig.dtype):
                # We need to infer the new categories.
                dtype = None
            else:
                dtype = self._orig.dtype
            result = Series(result, dtype=dtype, index=data.index, name=self._orig.name)
        return result


class ObjectProxy(PandasArray):
    _str = CachedAccessor("str", ObjectArrayMethods)
