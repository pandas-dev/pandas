import pytest

import pandas as pd
import pandas._testing as tm


@pytest.fixture
def input_dataframe():
    """Create the input dataframe"""

    # create input data
    input_dict = {
        "col1": [1, 2, 3, 4],
        "col2": ["a", "b", "c", "d"],
        "col3": [1.5, 2.5, 3.5, 4.5],
        "col4": ["cat1", "cat2", "cat3", "cat4"],
        "col5": ["obj1", "obj2", "obj3", "obj4"],
    }

    # explicitly cast columns as category and order them
    input_df = pd.DataFrame(data=input_dict).astype(
        {"col2": "category", "col4": "category"}
    )
    input_df["col2"] = input_df["col2"].cat.reorder_categories(
        ["a", "b", "c", "d"], ordered=True
    )
    input_df["col4"] = input_df["col4"].cat.reorder_categories(
        ["cat1", "cat2", "cat3", "cat4"], ordered=True
    )

    return input_df


@pytest.fixture
def expected_dataframe():
    """create the expected dataframe"""

    # create expected dataframe
    expected_dict = {
        "col1": [1, 2, 3, 4],
        "col2": ["a", "b", "c", "z"],
        "col3": [1.5, 2.5, 3.5, 4.5],
        "col4": ["cat1", "catX", "cat3", "cat4"],
        "col5": ["obj9", "obj2", "obj3", "obj4"],
    }

    # explicitly cast columns as category and order them
    expected_df = pd.DataFrame(data=expected_dict).astype(
        {"col2": "category", "col4": "category"}
    )
    expected_df["col2"] = expected_df["col2"].cat.reorder_categories(
        ["a", "b", "c", "z"], ordered=True
    )
    expected_df["col4"] = expected_df["col4"].cat.reorder_categories(
        ["cat1", "catX", "cat3", "cat4"], ordered=True
    )

    return expected_df


def test_replace_values_scalar(input_dataframe, expected_dataframe):

    # replace values in input dataframe
    input_df = input_dataframe.replace("d", "z")
    input_df = input_df.replace("obj1", "obj9")
    input_df = input_df.replace("cat2", "catX")

    tm.assert_frame_equal(input_df, expected_dataframe)
