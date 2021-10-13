import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


@pytest.fixture
def education_df():
    return pd.DataFrame(
        {
            "gender": ["male", "male", "female", "male", "female", "male"],
            "education": ["low", "medium", "high", "low", "high", "low"],
            "country": ["US", "FR", "US", "FR", "FR", "US"],
        }
    )


@pytest.mark.parametrize("normalize", [True, False])
@pytest.mark.parametrize(
    "sort, ascending",
    [
        (False, None),
        (True, True),
        (True, False),
    ],
)
def test_basic(education_df, normalize, sort, ascending):
    # gh43564
    # compare against apply with DataFrame value_counts
    def compute_proportions(df, var):
        return df[var].value_counts(normalize=normalize, sort=sort, ascending=ascending)

    expected = education_df.groupby("country").apply(
        compute_proportions, var=["gender", "education"]
    )
    result = education_df.groupby("country")[["gender", "education"]].value_counts(
        normalize=normalize, sort=sort, ascending=ascending
    )

    tm.assert_series_equal(result, expected)

    # compare against Series value_counts
    education_df["both"] = education_df["gender"] + "-" + education_df["education"]
    expected = education_df.groupby("country")["both"].value_counts(
        normalize=normalize, sort=sort, ascending=ascending
    )

    assert np.array_equal(result.values, expected.values)
