import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm

RESULT_NAME = "count"


def test_axis(animals_df):
    gp = animals_df.groupby([0, 0, 0], axis=1)
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


def test_basic(education_df):
    # gh43564
    result = education_df.groupby("country")[["gender", "education"]].value_counts(
        normalize=True
    )
    expected = pd.Series(
        name="count",
        data=[0.5, 0.25, 0.25, 0.5, 0.5],
        index=pd.MultiIndex.from_tuples(
            [
                ("FR", "male", "low"),
                ("FR", "female", "high"),
                ("FR", "male", "medium"),
                ("US", "female", "high"),
                ("US", "male", "low"),
            ],
            names=["country", "gender", "education"],
        ),
    )
    tm.assert_series_equal(result, expected)


def _frame_value_counts(df, keys, normalize, sort, ascending):
    return df[keys].value_counts(normalize=normalize, sort=sort, ascending=ascending)


@pytest.mark.parametrize("groupby", ["column", "array", "function"])
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
def test_against_frame_and_seriesgroupby(
    education_df, groupby, normalize, sort, ascending, as_index, frame
):
    # test all parameters:
    # - Use column, array or function as by= parameter
    # - Whether or not to normalize
    # - Whether or not to sort and how
    # - Whether or not to use the groupby as an index
    # - 3-way compare against:
    #   - apply with :meth:`~DataFrame.value_counts`
    #   - `~SeriesGroupBy.value_counts` (apart from certain cases where it crashes)
    by = {
        "column": "country",
        "array": education_df["country"].values,
        "function": lambda x: education_df["country"][x] == "US",
    }[groupby]

    gp = education_df.groupby(by=by, as_index=as_index)
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
            expected = expected.reset_index().rename({0: "count"}, axis=1)
            if groupby == "column":
                expected = expected.rename({"level_0": "country"}, axis=1)
                expected["country"] = np.where(expected["country"], "US", "FR")
            elif groupby == "function":
                expected["level_0"] = expected["level_0"] == 1
            else:
                expected["level_0"] = np.where(expected["level_0"], "US", "FR")
            tm.assert_frame_equal(result, expected)
    elif groupby == "column" or as_index:
        # (otherwise SeriesGroupby crashes)
        # compare against SeriesGroupBy value_counts
        education_df["both"] = education_df["gender"] + "-" + education_df["education"]
        expected = gp["both"].value_counts(
            normalize=normalize, sort=sort, ascending=ascending
        )
        if as_index:
            index_frame = expected.index.to_frame(index=False)
            index_frame["gender"] = index_frame["both"].str.split("-").str.get(0)
            index_frame["education"] = index_frame["both"].str.split("-").str.get(1)
            del index_frame["both"]
            index_frame = index_frame.rename({0: None}, axis=1)
            expected.index = pd.MultiIndex.from_frame(index_frame)
            expected.name = RESULT_NAME
            tm.assert_series_equal(result, expected)
        else:
            expected.insert(1, "gender", expected["both"].str.split("-").str.get(0))
            expected.insert(2, "education", expected["both"].str.split("-").str.get(1))
            del expected["both"]
            tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("normalize", [True, False])
@pytest.mark.parametrize(
    "sort, ascending, expected_rows, expected_count, expected_group_size",
    [
        (False, None, [0, 1, 2, 3, 4], [1, 1, 1, 2, 1], [1, 3, 1, 3, 1]),
        (True, False, [4, 3, 1, 2, 0], [1, 2, 1, 1, 1], [1, 3, 3, 1, 1]),
        (True, True, [4, 1, 3, 2, 0], [1, 1, 2, 1, 1], [1, 3, 3, 1, 1]),
    ],
)
def test_compound(
    education_df,
    normalize,
    sort,
    ascending,
    expected_rows,
    expected_count,
    expected_group_size,
):
    # Multiple groupby keys and as_index=False
    gp = education_df.groupby(["country", "gender"], as_index=False, sort=False)
    result = gp["education"].value_counts(
        normalize=normalize, sort=sort, ascending=ascending
    )
    expected = pd.DataFrame()
    for column in ["country", "gender", "education"]:
        expected[column] = [education_df[column][row] for row in expected_rows]
    expected["count"] = expected_count
    if normalize:
        expected["count"] /= expected_group_size
    tm.assert_frame_equal(result, expected)


@pytest.fixture
def animals_df():
    return pd.DataFrame(
        {"key": [1, 1, 1, 1], "num_legs": [2, 4, 4, 6], "num_wings": [2, 0, 0, 0]},
        index=["falcon", "dog", "cat", "ant"],
    )


@pytest.mark.parametrize(
    "sort, ascending, normalize, expected_data, expected_index",
    [
        (False, None, False, [1, 2, 1], [(1, 1, 1), (2, 4, 6), (2, 0, 0)]),
        (True, True, False, [1, 1, 2], [(1, 1, 1), (2, 6, 4), (2, 0, 0)]),
        (True, False, False, [2, 1, 1], [(1, 1, 1), (4, 2, 6), (0, 2, 0)]),
        (True, False, True, [0.5, 0.25, 0.25], [(1, 1, 1), (4, 2, 6), (0, 2, 0)]),
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
            expected_index, names=["key", "num_legs", "num_wings"]
        ),
    )
    tm.assert_series_equal(result_frame, expected)

    expected.name = RESULT_NAME
    result_frame_groupby = animals_df.groupby("key").value_counts(
        sort=sort, ascending=ascending, normalize=normalize
    )

    tm.assert_series_equal(result_frame_groupby, expected)


@pytest.fixture
def names_with_nulls_df(nulls_fixture):
    return pd.DataFrame(
        {
            "key": [1, 1, 1, 1],
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
                [(1, 1), ("Beth", "John"), ("Louise", "Smith")],
                names=["key", "first_name", "middle_name"],
            ),
        ),
        (
            False,
            [1, 1, 1, 1],
            pd.MultiIndex(
                levels=[
                    pd.Index([1]),
                    pd.Index(["Anne", "Beth", "John"]),
                    pd.Index(["Louise", "Smith", np.nan]),
                ],
                codes=[[0, 0, 0, 0], [0, 1, 2, 2], [2, 0, 1, 2]],
                names=["key", "first_name", "middle_name"],
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

    expected.name = RESULT_NAME
    result_frame_groupby = names_with_nulls_df.groupby("key").value_counts(
        dropna=dropna, normalize=normalize
    )

    tm.assert_series_equal(result_frame_groupby, expected)


@pytest.mark.parametrize("as_index", [False, True])
@pytest.mark.parametrize("observed", [False, True])
@pytest.mark.parametrize("normalize", [False, True])
def test_categorical(education_df, as_index, observed, normalize):
    # Test categorical data whether or not observed
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
