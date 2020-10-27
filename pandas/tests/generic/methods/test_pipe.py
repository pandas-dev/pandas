import pytest

from pandas import DataFrame, Series
import pandas._testing as tm


class TestPipe:
    @pytest.mark.parametrize("klass", [Series, DataFrame])
    def test_pipe(self, klass):
        obj = DataFrame({"A": [1, 2, 3]})
        expected = DataFrame({"A": [1, 4, 9]})
        if klass is Series:
            obj = obj["A"]
            expected = expected["A"]

        f = lambda x, y: x ** y
        result = obj.pipe(f, 2)
        tm.assert_equal(result, expected)

    @pytest.mark.parametrize("klass", [Series, DataFrame])
    def test_pipe_tuple(self, klass):
        obj = DataFrame({"A": [1, 2, 3]})
        if klass is Series:
            obj = obj["A"]

        f = lambda x, y: y
        result = obj.pipe((f, "y"), 0)
        tm.assert_equal(result, obj)

    @pytest.mark.parametrize("klass", [Series, DataFrame])
    def test_pipe_tuple_error(self, klass):
        obj = DataFrame({"A": [1, 2, 3]})
        if klass is Series:
            obj = obj["A"]

        f = lambda x, y: y
        with pytest.raises(ValueError):
            obj.pipe((f, "y"), x=1, y=0)
