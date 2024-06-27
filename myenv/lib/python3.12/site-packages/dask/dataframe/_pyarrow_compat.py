from __future__ import annotations

import copyreg
import warnings

import pandas as pd
from packaging.version import Version

try:
    import pyarrow as pa

    pa_version = pa.__version__
    if Version(pa_version) < Version("14.0.1"):
        try:
            import pyarrow_hotfix  # noqa: F401

            warnings.warn(
                "Minimal version of pyarrow will soon be increased to 14.0.1. "
                f"You are using {pa_version}. Please consider upgrading.",
                FutureWarning,
            )
        except ImportError:
            warnings.warn(
                f"You are using pyarrow version {pa_version} which is "
                "known to be insecure. See https://www.cve.org/CVERecord?id=CVE-2023-47248 "
                "for further details. Please upgrade to pyarrow>=14.0.1 "
                "or install pyarrow-hotfix to patch your current version."
            )
except ImportError:
    pa = None

from dask.dataframe._compat import PANDAS_GE_150, PANDAS_GE_200

# Pickling of pyarrow arrays is effectively broken - pickling a slice of an
# array ends up pickling the entire backing array.
#
# See https://issues.apache.org/jira/browse/ARROW-10739
#
# This comes up when using pandas `string[pyarrow]` dtypes, which are backed by
# a `pyarrow.StringArray`.  To fix this, we register a *global* override for
# pickling `ArrowStringArray` or `ArrowExtensionArray` types (where available).
# We do this at the pandas level rather than the pyarrow level for efficiency reasons
# (a pandas ArrowStringArray may contain many small pyarrow StringArray objects).
#
# The implementation here is based on https://github.com/pandas-dev/pandas/pull/49078
# which is included in pandas=2+. We can remove all this once Dask's minimum
# supported pandas version is at least 2.0.0.


def rebuild_arrowextensionarray(type_, chunks):
    array = pa.chunked_array(chunks)
    return type_(array)


def reduce_arrowextensionarray(x):
    return (rebuild_arrowextensionarray, (type(x), x._data.combine_chunks()))


# `pandas=2` includes efficient serialization of `pyarrow`-backed extension arrays.
# See https://github.com/pandas-dev/pandas/pull/49078 for details.
# We only need to backport efficient serialization for `pandas<2`.
if pa is not None and not PANDAS_GE_200:
    if PANDAS_GE_150:
        # Applies to all `pyarrow`-backed extension arrays (e.g. `string[pyarrow]`, `int64[pyarrow]`)
        for type_ in [pd.arrays.ArrowExtensionArray, pd.arrays.ArrowStringArray]:
            copyreg.dispatch_table[type_] = reduce_arrowextensionarray
    else:
        # Only `string[pyarrow]` is implemented, so just patch that
        copyreg.dispatch_table[pd.arrays.ArrowStringArray] = reduce_arrowextensionarray
