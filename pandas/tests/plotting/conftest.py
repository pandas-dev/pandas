import numpy as np
import pytest

from pandas import (
    DataFrame,
    to_datetime,
)


@pytest.fixture
def hist_df():
    n = 100
    np_random = np.random.RandomState(42)
    gender = np_random.choice(["Male", "Female"], size=n)
    classroom = np_random.choice(["A", "B", "C"], size=n)

    hist_df = DataFrame(
        {
            "gender": gender,
            "classroom": classroom,
            "height": np.random.normal(66, 4, size=n),
            "weight": np.random.normal(161, 32, size=n),
            "category": np.random.randint(4, size=n),
            "datetime": to_datetime(
                np.random.randint(
                    812419200000000000,
                    819331200000000000,
                    size=n,
                    dtype=np.int64,
                )
            ),
        }
    )
    return hist_df
