import operator

import pytest

from pandas import Series


@pytest.fixture(name="dtype")
def fixture_dtype():
    """A fixture providing the ExtensionDtype to validate."""
    raise NotImplementedError


@pytest.fixture(name="data")
def fixture_data():
    """
    Length-100 array for this type.

    * data[0] and data[1] should both be non missing
    * data[0] and data[1] should not be equal
    """
    raise NotImplementedError


@pytest.fixture(name="data_for_twos")
def fixture_data_for_twos():
    """Length-100 array in which all the elements are two."""
    raise NotImplementedError


@pytest.fixture(name="data_missing")
def fixture_data_missing():
    """Length-2 array with [NA, Valid]"""
    raise NotImplementedError


@pytest.fixture(name="all_data", params=["data", "data_missing"])
def fixture_all_data(request, data, data_missing):
    """Parametrized fixture giving 'data' and 'data_missing'"""
    if request.param == "data":
        return data
    elif request.param == "data_missing":
        return data_missing


@pytest.fixture(name="data_repeated")
def fixture_data_repeated(data):
    """
    Generate many datasets.

    Parameters
    ----------
    data : fixture implementing `data`

    Returns
    -------
    Callable[[int], Generator]:
        A callable that takes a `count` argument and
        returns a generator yielding `count` datasets.
    """

    def gen(count):
        for _ in range(count):
            yield data

    return gen


@pytest.fixture(name="data_for_sorting")
def fixture_data_for_sorting():
    """
    Length-3 array with a known sort order.

    This should be three items [B, C, A] with
    A < B < C
    """
    raise NotImplementedError


@pytest.fixture(name="data_missing_for_sorting")
def fixture_data_missing_for_sorting():
    """
    Length-3 array with a known sort order.

    This should be three items [B, NA, A] with
    A < B and NA missing.
    """
    raise NotImplementedError


@pytest.fixture(name="na_cmp")
def fixture_na_cmp():
    """
    Binary operator for comparing NA values.

    Should return a function of two arguments that returns
    True if both arguments are (scalar) NA for your type.

    By default, uses ``operator.is_``
    """
    return operator.is_


@pytest.fixture(name="na_value")
def fixture_na_value():
    """The scalar missing value for this type. Default 'None'"""
    return None


@pytest.fixture(name="data_for_grouping")
def fixture_data_for_grouping():
    """
    Data for factorization, grouping, and unique tests.

    Expected to be like [B, B, NA, NA, A, A, B, C]

    Where A < B < C and NA is missing
    """
    raise NotImplementedError


@pytest.fixture(name="box_in_series", params=[True, False])
def fixture_box_in_series(request):
    """Whether to box the data in a Series"""
    return request.param


@pytest.fixture(
    name="groupby_apply_op",
    params=[
        lambda x: 1,
        lambda x: [1] * len(x),
        lambda x: Series([1] * len(x)),
        lambda x: x,
    ],
    ids=["scalar", "list", "series", "object"],
)
def fixture_groupby_apply_op(request):
    """
    Functions to test groupby.apply().
    """
    return request.param


@pytest.fixture(name="as_frame", params=[True, False])
def fixture_as_frame(request):
    """
    Boolean fixture to support Series and Series.to_frame() comparison testing.
    """
    return request.param


@pytest.fixture(name="as_series", params=[True, False])
def fixture_as_series(request):
    """
    Boolean fixture to support arr and Series(arr) comparison testing.
    """
    return request.param


@pytest.fixture(name="use_numpy", params=[True, False])
def fixture_use_numpy(request):
    """
    Boolean fixture to support comparison testing of ExtensionDtype array
    and numpy array.
    """
    return request.param


@pytest.fixture(name="fillna_method", params=["ffill", "bfill"])
def fixture_fillna_method(request):
    """
    Parametrized fixture giving method parameters 'ffill' and 'bfill' for
    Series.fillna(method=<method>) testing.
    """
    return request.param


@pytest.fixture(name="as_array", params=[True, False])
def fixture_as_array(request):
    """
    Boolean fixture to support ExtensionDtype _from_sequence method testing.
    """
    return request.param


@pytest.fixture(name="invalid_scalar")
def fixture_invalid_scalar(data):
    """
    A scalar that *cannot* be held by this ExtensionArray.

    The default should work for most subclasses, but is not guaranteed.

    If the array can hold any item (i.e. object dtype), then use pytest.skip.
    """
    return object.__new__(object)
