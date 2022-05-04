import pytest


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

    if func == "_set_name":
        kwargs["name"] = "hello"

    with pytest.raises(ValueError, match=msg):
        getattr(string_series, func)(**kwargs)

@pytest.mark.parametrize('keyword', ('copy', 'fastpath', 'takeable', 'clear', 'verify_is_copy', 'inplace', 'allow_duplicates', 'index', 'as_index', 'sort', 'observed', 'dropna', 'ignore_index', 'verify_integrity', 'keep_shape', 'keep_equal', 'inplace', 'sort_remaining' , 'convert_dtype', 'show_counts', 'deep', 'infer_objects', 'convert_string', 'convert_integer', 'convert_boolean', 'convert_floating', 'normalize'))
def test_set_index_validation(string_series, func, keyword):
    msg = 'For argument "{}" expected type bool'.format(keyword)
    kwargs = {keyword: 'hello'}
    with pytest.raises(ValueError, match=msg):
        getattr(string_series, func)(**kwargs)