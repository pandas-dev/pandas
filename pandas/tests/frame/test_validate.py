import pytest

from pandas.core.frame import DataFrame


@pytest.fixture
def dataframe():
    return DataFrame({"a": [1, 2], "b": [3, 4]})


class TestDataFrameValidate:
    """Tests for error handling related to data types of method arguments."""

    @pytest.mark.parametrize(
        "func",
        [
            "query",
            "eval",
            "set_index",
            "reset_index",
            "dropna",
            "drop_duplicates",
            "sort_values",
        ],
    )
    @pytest.mark.parametrize("inplace", [1, "True", [1, 2, 3], 5.0])
    def test_validate_bool_args(self, dataframe, func, inplace):
        msg = 'For argument "inplace" expected type bool'
        kwargs = {"inplace": inplace}

        if func == "query":
            kwargs["expr"] = "a > b"
        elif func == "eval":
            kwargs["expr"] = "a + b"
        elif func == "set_index":
            kwargs["keys"] = ["a"]
        elif func == "sort_values":
            kwargs["by"] = ["a"]

        with pytest.raises(ValueError, match=msg):
            getattr(dataframe, func)(**kwargs)

@pytest.mark.parametrize('keyword', ('nan_as_null', 'allow_copy', 'ignore_width', 'index', 'index_names', 'show_dimensions', 'copy', 'inplace', 'reauth', 'auth_local_webserver', 'progress_bar', 'verify_integrity', 'write_index', 'bold_rows', 'escape', 'notebook', 'render_links', 'deep', 'takeable', 'drop', 'append', 'ignore_index', 'sort_remaining', 'normalize', 'ascending', 'dropna', 'keep_shape', 'keep_equal', 'overwrite', 'as_index', 'observed', 'sort', 'raw', 'left_index', 'right_index', 'numeric_only', 'skipna'))
def test_set_index_validation(dataframe, func, keyword):
    msg = 'For argument "{}" expected type bool'.format(keyword)
    kwargs = {keyword: 'hello'}
    with pytest.raises(ValueError, match=msg):
        getattr(dataframe, func)(**kwargs)