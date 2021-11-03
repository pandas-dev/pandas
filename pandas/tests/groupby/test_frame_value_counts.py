import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm

RESULT_NAME = "count"


def test_axis(animals_df):
    gp = animals_df.groupby([0, 0], axis=1)
    with pytest.raises(NotImplementedError, match="axis"):
        gp.value_counts()


@pytest.fixture
def education_df():
    return pd.DataFrame(
        {
            "gender": ["male", "male", "female", "male", "female", "male"],
            "education": ["low", "medium", "high", "low", "high", "low"],
            "country": ["US", "FR", "US", "FR", "FR", "FR"],
        }
    )


@pytest.fixture
def country_index():
    return ["U", "F", "U", "F", "F", "F"]


def _frame_value_counts(df, keys, normalize, sort, ascending):
    return df[keys].value_counts(normalize=normalize, sort=sort, ascending=ascending)


@pytest.mark.parametrize("column", [True, False])
@pytest.mark.parametrize("normalize", [True, False])
@pytest.mark.parametrize(
    "sort, ascending",
    [
        (False, None),
        (True, True),
        (True, False),
    ],
)
@pytest.mark.parametrize("as_index", [True, False])
@pytest.mark.parametrize("frame", [True, False])
def test_basic(
    education_df, country_index, column, normalize, sort, ascending, as_index, frame
):
    # gh43564 with added:
    # - Use column or index
    # - Whether or not to normalize
    # - Whether or not to sort and how
    # - Whether or not to use the groupby as an index
    # - 3-way compare against :meth:`~DataFrame.value_counts`
    #   and `~SeriesGroupBy.value_counts`
    gp = education_df.groupby("country" if column else country_index, as_index=as_index)
    result = gp[["gender", "education"]].value_counts(
        normalize=normalize, sort=sort, ascending=ascending
    )

    if frame:
        # compare against apply with DataFrame value_counts
        expected = gp.apply(
            _frame_value_counts, ["gender", "education"], normalize, sort, ascending
        )

        expected.name = RESULT_NAME
        if as_index:
            tm.assert_series_equal(result, expected)
        else:
            tm.assert_numpy_array_equal(result[RESULT_NAME].values, expected.values)
    elif column or as_index:
        # (otherwise SeriesGroupby crashes)
        # compare against SeriesGroupBy value_counts
        education_df["both"] = education_df["gender"] + "-" + education_df["education"]
        expected = gp["both"].value_counts(
            normalize=normalize, sort=sort, ascending=ascending
        )
        if as_index:
            tm.assert_numpy_array_equal(result.values, expected.values)
        else:
            tm.assert_numpy_array_equal(
                result[RESULT_NAME].values, expected[RESULT_NAME].values
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
def test_compound(education_df, normalize, sort, ascending):
    # Multiple groupby keys and as_index=False
    gp = education_df.groupby(["country", "gender"], as_index=False, sort=False)
    result = gp["education"].value_counts(
        normalize=normalize, sort=sort, ascending=ascending
    )

    if sort:
        # compare against apply with DataFrame value_counts
        expected_values = gp.apply(
            _frame_value_counts, "education", normalize, sort, ascending
        ).values
    elif normalize:
        expected_values = np.array([1.0, 1.0 / 3, 1.0, 2.0 / 3, 1.0])
    else:
        expected_values = np.array([1, 1, 1, 2, 1], dtype=np.int64)

    tm.assert_numpy_array_equal(result[RESULT_NAME].values, expected_values)


@pytest.fixture
def animals_df():
    return pd.DataFrame(
        {"num_legs": [2, 4, 4, 6], "num_wings": [2, 0, 0, 0]},
        index=["falcon", "dog", "cat", "ant"],
    )


@pytest.mark.parametrize(
    "sort, ascending, normalize, expected_data, expected_index",
    [
        (False, None, False, [1, 2, 1], [(2, 4, 6), (2, 0, 0)]),
        (True, True, False, [1, 1, 2], [(2, 6, 4), (2, 0, 0)]),
        (True, False, False, [2, 1, 1], [(4, 2, 6), (0, 2, 0)]),
        (True, False, False, [2, 1, 1], [(4, 2, 6), (0, 2, 0)]),
        (True, False, True, [0.5, 0.25, 0.25], [(4, 2, 6), (0, 2, 0)]),
    ],
)
def test_data_frame_value_counts(
    animals_df, sort, ascending, normalize, expected_data, expected_index
):
    # 3-way compare with :meth:`~DataFrame.value_counts`
    # Tests from frame/methods/test_value_counts.py
    result_frame = animals_df.value_counts(
        sort=sort, ascending=ascending, normalize=normalize
    )
    expected = pd.Series(
        data=expected_data,
        index=pd.MultiIndex.from_arrays(
            expected_index, names=["num_legs", "num_wings"]
        ),
    )
    tm.assert_series_equal(result_frame, expected)

    animals_df["key"] = 1

    result_frame_groupby = animals_df.groupby("key").value_counts(
        sort=sort, ascending=ascending, normalize=normalize
    )
    result_frame_groupby.reset_index(drop=True, level="key", inplace=True)
    result_frame_groupby.name = None
    tm.assert_series_equal(result_frame_groupby, expected)


@pytest.fixture
def names_with_nulls_df(nulls_fixture):
    return pd.DataFrame(
        {
            "first_name": ["John", "Anne", "John", "Beth"],
            "middle_name": ["Smith", nulls_fixture, nulls_fixture, "Louise"],
        },
    )


@pytest.mark.parametrize(
    "dropna, expected_data, expected_index",
    [
        (
            True,
            [1, 1],
            pd.MultiIndex.from_arrays(
                [("Beth", "John"), ("Louise", "Smith")],
                names=["first_name", "middle_name"],
            ),
        ),
        (
            False,
            [1, 1, 1, 1],
            pd.MultiIndex(
                levels=[
                    pd.Index(["Anne", "Beth", "John"]),
                    pd.Index(["Louise", "Smith", np.nan]),
                ],
                codes=[[0, 1, 2, 2], [2, 0, 1, 2]],
                names=["first_name", "middle_name"],
            ),
        ),
    ],
)
@pytest.mark.parametrize("normalize", [False, True])
def test_data_frame_value_counts_dropna(
    names_with_nulls_df, dropna, normalize, expected_data, expected_index
):
    # GH 41334
    # 3-way compare with :meth:`~DataFrame.value_counts`
    # Tests with nulls from frame/methods/test_value_counts.py
    result_frame = names_with_nulls_df.value_counts(dropna=dropna, normalize=normalize)
    expected = pd.Series(
        data=expected_data,
        index=expected_index,
    )
    if normalize:
        expected /= float(len(expected_data))

    tm.assert_series_equal(result_frame, expected)

    names_with_nulls_df["key"] = 1
    result_frame_groupby = names_with_nulls_df.groupby("key").value_counts(
        dropna=dropna, normalize=normalize
    )
    result_frame_groupby.reset_index(drop=True, level="key", inplace=True)
    result_frame_groupby.name = None

    tm.assert_series_equal(result_frame_groupby, expected)


@pytest.mark.parametrize("as_index", [False, True])
@pytest.mark.parametrize("observed", [False, True])
@pytest.mark.parametrize("normalize", [False, True])
def test_categorical(education_df, as_index, observed, normalize):
    gp = education_df.astype("category").groupby(
        "country", as_index=as_index, observed=observed
    )
    result = gp.value_counts(normalize=normalize)

    if observed:
        gp = education_df.groupby("country", as_index=as_index)
        expected = gp.value_counts(normalize=normalize)
        if as_index:
            tm.assert_numpy_array_equal(result.values, expected.values)
        else:
            tm.assert_numpy_array_equal(
                result[RESULT_NAME].values, expected[RESULT_NAME].values
            )
    else:
        if normalize:
            expected_values = np.array(
                [0.5, 0.25, 0.25, 0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0]
            )
        else:
            expected_values = np.array(
                [2, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0], dtype=np.int64
            )

        if as_index:
            tm.assert_numpy_array_equal(result.values, expected_values)
        else:
            tm.assert_numpy_array_equal(result[RESULT_NAME].values, expected_values)
