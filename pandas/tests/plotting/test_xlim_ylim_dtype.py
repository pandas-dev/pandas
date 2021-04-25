import pytest
import pandas as pd
from datetime import datetime as dt
import numpy as np

"""targeting #GH40781"""

# arbitrary df to test
df = pd.DataFrame(
    {
        "A" : [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        "B": [dt.now() for i in range(0, 10)]
    }
)


def test_axis_lim_invalid():
    global df
    """
    supplying ranges with str, datetime and
    mixed int,str and float,str dtype to raise
    ValueError. Valid dtypes are np.datetime64
    int, and float.
    """
    lim = [
        ["1", "3"],
        [dt.now(), dt.now()],
        [1, "2"],
        [0.1, "0.2"],
        [np.datetime64(dt.now()), dt.now()]
    ]
    for elem in lim:
        with pytest.raises(ValueError, match="`xlim` contains values"):
            df.plot.line(x="A", xlim=elem)


def test_axis_lim_valid():
    global df
    """
    supplying ranges with
    valid dtypes: np.datetime64
    int, and float, this test
    should not raise errors.
    """
    lim = [
        [1, 3],
        [np.datetime64(dt.now()), np.datetime64(dt.now())],
        [0.1, 0.2],
    ]
    for elem in lim:
        df.plot.line(x="A", ylim=elem)
