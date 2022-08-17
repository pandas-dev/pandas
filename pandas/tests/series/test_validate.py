import pytest

import pandas._testing as tm


@pytest.mark.parametrize(
    "func",
    [
        "reset_index",
        "_set_name",
        "sort_values",
        "sort_index",
        "rename",
        "dropna",
        "drop_duplicates",
    ],
)
@pytest.mark.parametrize("inplace", [1, "True", [1, 2, 3], 5.0])
def test_validate_bool_args(string_series, func, inplace):
    """Tests for error handling related to data types of method arguments."""
    msg = 'For argument "inplace" expected type bool'
    kwargs = {"inplace": inplace}

    warn_msg = "'inplace' keyword in Series.rename"
    warn = None
    if func == "_set_name":
        kwargs["name"] = "hello"
        warn = FutureWarning
    elif func == "rename":
        warn = FutureWarning

    with pytest.raises(ValueError, match=msg):
        with tm.assert_produces_warning(warn, match=warn_msg):
            getattr(string_series, func)(**kwargs)
