import numpy as np
import pytest

from pandas.errors import Pandas4Warning

from pandas import Index
import pandas._testing as tm


def test_invert_object_deprecated():
    # GH#51567 ~ on object dtype dispatches to the objects' own __invert__
    idx = Index([1, 2], dtype=object)

    msg = "__invert__ .* on object dtype is deprecated"
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        result = ~idx
    tm.assert_index_equal(result, Index([-2, -3], dtype=object))


def test_invert_object_bool_mask_deprecated():
    # GH#16873, GH#31035 a boolean mask cast to object silently becomes
    #  integers instead of being negated
    idx = Index(np.array([True, False], dtype=object))

    msg = "__invert__ .* on object dtype is deprecated"
    # Python itself deprecated ~ on bool in 3.12, so on new enough versions
    #  there is a second warning here that we don't care about.
    with tm.assert_produces_warning(
        Pandas4Warning, match=msg, raise_on_extra_warnings=False
    ):
        result = ~idx
    tm.assert_index_equal(result, Index(np.array([-2, -1], dtype=object)))


def test_invert_object_raising_still_warns():
    # GH#51567 elements that have no __invert__ warn before raising
    idx = Index([1.0, 2.0], dtype=object)

    msg = "__invert__ .* on object dtype is deprecated"
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        with pytest.raises(TypeError, match="bad operand type for unary ~"):
            ~idx
