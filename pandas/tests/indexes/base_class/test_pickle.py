from pandas import Index
import pandas._testing as tm


def test_pickle_preserves_object_dtype(tmp_path):
    # GH#43188, GH#43155 don't infer numeric dtype
    index = Index([1, 2, 3], dtype=object)

    result = tm.round_trip_pickle(index, tmp_path)
    assert result.dtype == object
    tm.assert_index_equal(index, result)
