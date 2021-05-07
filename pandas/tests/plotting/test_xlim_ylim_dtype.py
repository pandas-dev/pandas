from datetime import datetime as dt

import numpy as np
import pytest

pytest.importorskip("matplotlib")
from pandas import DataFrame

pytestmark = pytest.mark.slow


@pytest.mark.parametrize("df, lim", [
    (
        DataFrame(
            {
                "A" : [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                "B": [dt.now() for i in range(0, 10)]
            }
        ),
        [
            ["1", "3"],
            [dt.now(), dt.now()],
            [1, "2"],
            [0.1, "0.2"],
            [np.datetime64(dt.now()), dt.now()]
        ]
    )
]
)
def test_axis_lim_invalid(df, lim):
    for elem in lim:
        with pytest.raises(ValueError, match="`xlim` contains values"):
            df.plot.line(x="A", xlim=elem)


@pytest.mark.parametrize("df, lim", [
    (
        DataFrame(
            {
                "A" : [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                "B": [dt.now() for i in range(0, 10)]
            }
        ),
        [
            [1, 3],
            [0.1, 0.2],
            [np.datetime64(dt.now()), np.datetime64(dt.now())]
        ]
    )
]
)
def test_axis_lim_valid(df, lim):
    for elem in lim:
        df.plot.line(x="A", xlim=elem)
