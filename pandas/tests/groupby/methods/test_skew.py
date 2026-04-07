import numpy as np
import pytest

from pandas.errors import Pandas4Warning

import pandas as pd
import pandas._testing as tm


def test_groupby_skew_equivalence():
    # Test that the groupby skew method (which uses libgroupby.group_skew)
    #  matches the results of operating group-by-group (which uses nanops.nanskew)
    nrows = 1000
    ngroups = 3
    ncols = 2
    nan_frac = 0.05

    arr = np.random.default_rng(2).standard_normal((nrows, ncols))
    arr[np.random.default_rng(2).random(nrows) < nan_frac] = np.nan

    df = pd.DataFrame(arr)
    grps = np.random.default_rng(2).integers(0, ngroups, size=nrows)
    gb = df.groupby(grps)

    result = gb.skew()

    grpwise = [grp.skew().to_frame(i).T for i, grp in gb]
    expected = pd.concat(grpwise, axis=0)
    expected.index = expected.index.astype(result.index.dtype)  # 32bit builds
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("klass", ["SeriesGroupBy", "DataFrameGroupBy"])
def test_skew_kwargs_deprecated(klass):
    # GH#50407
    df = pd.DataFrame({"a": [1, 1, 2], "b": [1.0, 2.0, 3.0]})
    if klass == "SeriesGroupBy":
        gb = df.groupby("a")["b"]
        msg = "Passing additional arguments to SeriesGroupBy.skew"
    else:
        gb = df.groupby("a")
        msg = "Passing additional arguments to DataFrameGroupBy.skew"

    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        with pytest.raises(TypeError, match="got an unexpected keyword argument"):
            gb.skew(foo=1)
