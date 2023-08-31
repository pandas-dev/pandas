import pandas as pd
import pandas._testing as tm
from pandas.api.types import is_list_like
from pandas.core.groupby.base import transformation_kernels

# There is no Series.cumcount or DataFrame.cumcount
series_transform_kernels = [
    x for x in sorted(transformation_kernels) if x != "cumcount"
]
frame_transform_kernels = [x for x in sorted(transformation_kernels) if x != "cumcount"]


def transform_obj(obj, func, *args, axis=0, series_ops_only=False, **kwargs):
    """helper function to ease use of series_ops_only and deprecation warning."""
    if series_ops_only:
        result = obj.transform(
            func, axis, *args, series_ops_only=series_ops_only, **kwargs
        )
    elif isinstance(obj, pd.DataFrame) and not is_list_like(func):
        result = obj.transform(func, axis, *args, **kwargs)
    elif isinstance(func, str):
        result = obj.transform(func, axis, *args, **kwargs)
    else:
        cls_name = type(obj).__name__
        msg = (
            f"{cls_name}.transform will in the future only operate on "
            "whole series. Set series_ops_only = True to opt into the new behavior "
            f"or use {cls_name}.map to continue operating on series elements."
        )
        with tm.assert_produces_warning(FutureWarning, match=msg):
            result = obj.transform(
                func, axis, *args, series_ops_only=series_ops_only, **kwargs
            )
    return result
