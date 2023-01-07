import numpy as np
import pytest  # noqa

import pandas as pd
import pandas._testing as tm


class TestCaseWhen:
    def test_case_when_multiple_conditions_callable(self):
        df = pd.DataFrame(dict(a=[1, 2, 3], b=[4, 5, 6]))
        result = df.assign(
            new_column=pd.case_when(
                lambda x: x.a == 1,
                1,
                lambda x: (x.a > 1) & (x.b == 5),
                2,
            )
        )
        expected = df.assign(new_column=[1, 2, np.nan])
        tm.assert_frame_equal(result, expected)

    def test_case_when_multiple_conditions_array_series(self):
        df = pd.DataFrame(dict(a=[1, 2, 3], b=[4, 5, 6]))
        result = df.assign(
            new_column=pd.case_when(
                [True, False, False],
                1,
                pd.Series([False, True, False]),
                2,
            )
        )
        expected = df.assign(new_column=[1, 2, np.nan])
        tm.assert_frame_equal(result, expected)

    def test_case_when_multiple_conditions_callable_default(self):
        df = pd.DataFrame(dict(a=[1, 2, 3], b=[4, 5, 6]))
        result = df.assign(
            new_column=pd.case_when(
                lambda x: x.a == 1,
                1,
                lambda x: (x.a > 1) & (x.b == 5),
                2,
                default=-1,
            )
        )
        expected = df.assign(new_column=[1, 2, -1])
        tm.assert_frame_equal(result, expected)

    def test_case_when_multiple_conditions_callable_default_series(self):
        df = pd.DataFrame(dict(a=[1, 2, 3], b=[4, 5, 6]))
        result = df.assign(
            new_column=pd.case_when(
                lambda x: x.a == 1,
                1,
                lambda x: (x.a > 1) & (x.b == 5),
                2,
                default=df.b,
            )
        )
        expected = df.assign(new_column=[1, 2, 6])
        tm.assert_frame_equal(result, expected)
