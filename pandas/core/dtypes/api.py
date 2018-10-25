# flake8: noqa

import sys

from .common import (  # categorical; interval; datetimelike; string-like; sparse; numeric types; like
    is_array_like, is_bool, is_bool_dtype, is_categorical,
    is_categorical_dtype, is_complex, is_complex_dtype,
    is_datetime64_any_dtype, is_datetime64_dtype, is_datetime64_ns_dtype,
    is_datetime64tz_dtype, is_datetimetz, is_dict_like, is_dtype_equal,
    is_extension_type, is_file_like, is_float, is_float_dtype, is_hashable,
    is_int64_dtype, is_integer, is_integer_dtype, is_interval,
    is_interval_dtype, is_iterator, is_list_like, is_named_tuple, is_number,
    is_numeric_dtype, is_object_dtype, is_period, is_period_dtype, is_re,
    is_re_compilable, is_scalar, is_signed_integer_dtype, is_sparse,
    is_string_dtype, is_timedelta64_dtype, is_timedelta64_ns_dtype,
    is_unsigned_integer_dtype, pandas_dtype
)

# deprecated
m = sys.modules['pandas.core.dtypes.api']

for t in ['is_any_int_dtype', 'is_floating_dtype', 'is_sequence']:

    def outer(t=t):

        def wrapper(arr_or_dtype):
            import warnings
            import pandas
            warnings.warn("{t} is deprecated and will be "
                          "removed in a future version".format(t=t),
                          FutureWarning, stacklevel=3)
            return getattr(pandas.core.dtypes.common, t)(arr_or_dtype)
        return wrapper

    setattr(m, t, outer(t))

del sys, m, t, outer
