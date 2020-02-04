import pytest

from pandas import Timestamp


@pytest.mark.parametrize("kwargs", [{}, {"year": 2020}, {"year": 2020, "month": 1}])
def test_timestamp_constructor_missing_keyword(kwargs):
    # GH 31200
    msg = r"function missing required argument .* \(pos .\)"
    with pytest.raises(TypeError, match=msg):
        Timestamp(**kwargs)
