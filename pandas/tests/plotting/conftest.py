import numpy as np
import pytest

import pandas as pd


@pytest.fixture(autouse=True)
def autouse_mpl_cleanup(mpl_cleanup):
    pass


@pytest.fixture(autouse=True)
def restore_plot_params():
    # plot_params is process-global, so a test that sets it changes how every
    # later test in the same process plots datetime indexes
    original = dict(pd.plotting.plot_params)
    yield
    pd.plotting.plot_params.clear()
    pd.plotting.plot_params.update(original)


@pytest.fixture
def hist_df():
    n = 50
    rng = np.random.default_rng(10)
    gender = rng.choice(["Male", "Female"], size=n)
    classroom = rng.choice(["A", "B", "C"], size=n)

    hist_df = pd.DataFrame(
        {
            "gender": gender,
            "classroom": classroom,
            "height": rng.normal(66, 4, size=n),
            "weight": rng.normal(161, 32, size=n),
            "category": rng.integers(4, size=n),
            "datetime": pd.to_datetime(
                rng.integers(
                    812419200000000000,
                    819331200000000000,
                    size=n,
                    dtype=np.int64,
                )
            ),
        }
    )
    return hist_df
