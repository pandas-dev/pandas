from pandas.core.dtypes.cast import construct_1d_arraylike_from_scalar
from pandas.core.dtypes.dtypes import CategoricalDtype

from pandas import Categorical
from pandas.util import testing as tm


def test_cast_1d_array_like_from_scalar_categorical():
    # see gh-19565
    #
    # Categorical result from scalar did not maintain
    # categories and ordering of the passed dtype.
    cats = ["a", "b", "c"]
    cat_type = CategoricalDtype(categories=cats, ordered=False)
    expected = Categorical(["a", "a"], categories=cats)

    result = construct_1d_arraylike_from_scalar("a", len(expected), cat_type)
    tm.assert_categorical_equal(
        result, expected, check_category_order=True, check_dtype=True
    )
