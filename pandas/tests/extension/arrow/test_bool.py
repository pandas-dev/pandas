import numpy as np
import pytest

import pandas as pd
from pandas.tests.extension import base
import pandas.util.testing as tm

pytest.importorskip('pyarrow', minversion="0.10.0")

from .bool import ArrowBoolArray, ArrowBoolDtype  # isort:skip


@pytest.fixture
def dtype():
    return ArrowBoolDtype()


@pytest.fixture
def data():
    return ArrowBoolArray.from_scalars(np.random.randint(0, 2, size=100,
                                                         dtype=bool))


@pytest.fixture
def data_missing():
    return ArrowBoolArray.from_scalars([None, True])


class BaseArrowTests:
    pass


class TestDtype(BaseArrowTests, base.BaseDtypeTests):
    def test_array_type_with_arg(self, data, dtype):
        pytest.skip("GH-22666")


class TestInterface(BaseArrowTests, base.BaseInterfaceTests):
    def test_repr(self, data):
        raise pytest.skip("TODO")


class TestConstructors(BaseArrowTests, base.BaseConstructorsTests):
    def test_from_dtype(self, data):
        pytest.skip("GH-22666")

    # seems like some bug in isna on empty BoolArray returning floats.
    @pytest.mark.xfail(reason='bad is-na for empty data')
    def test_from_sequence_from_cls(self, data):
        super().test_from_sequence_from_cls(data)


class TestReduce(base.BaseNoReduceTests):
    def test_reduce_series_boolean(self):
        pass


class TestReduceBoolean(base.BaseBooleanReduceTests):
    pass


def test_is_bool_dtype(data):
    assert pd.api.types.is_bool_dtype(data)
    assert pd.core.common.is_bool_indexer(data)
    s = pd.Series(range(len(data)))
    result = s[data]
    expected = s[np.asarray(data)]
    tm.assert_series_equal(result, expected)


class TestReshape:

    @pytest.mark.parametrize('shape', [(1, -1), (-1, 1)])
    def test_ravel(self, data, shape):
        dim2 = data.reshape(shape)
        assert dim2.ndim == 2

        assert dim2.ravel().shape == data.shape

    def test_transpose(self, data):
        assert data.T.shape == data.shape

        rowlike = data.reshape(-1, 1)
        collike = data.reshape(1, -1)

        rt1 = rowlike.T
        rt2 = rowlike.transpose()

        ct1 = collike.T
        ct2 = collike.transpose()

        assert ct1.shape == ct2.shape == rowlike.shape
        assert rt1.shape == rt2.shape == collike.shape

    def test_swapaxes(self, data):
        rowlike = data.reshape(-1, 1)
        collike = data.reshape(1, -1)

        for arr in [data, rowlike, collike]:
            assert arr.swapaxes(0, 0).shape == arr.shape

        assert rowlike.swapaxes(0, 1).shape == collike.shape
        assert rowlike.swapaxes(1, 0).shape == collike.shape
        assert collike.swapaxes(0, 1).shape == rowlike.shape
        assert collike.swapaxes(1, 0).shape == rowlike.shape

    def test_reshape(self, data):
        size = len(data)

        collike1 = data.reshape(1, -1)
        collike2 = data.reshape((1, -1))
        collike3 = data.reshape(1, size)
        collike4 = data.reshape((1, size))
        collike5 = collike1.reshape(collike1.shape)
        for collike in [collike1, collike2, collike3, collike4, collike5]:
            assert collike.shape == (1, size)
            assert collike.ndim == 2

        rowlike1 = data.reshape(-1, 1)
        rowlike2 = data.reshape((-1, 1))
        rowlike3 = data.reshape(size, 1)
        rowlike4 = data.reshape((size, 1))
        rowlike5 = rowlike1.reshape(rowlike1.shape)
        rowlike6 = collike1.reshape(-1, 1)
        for collike in [rowlike1, rowlike2, rowlike3, rowlike4,
                        rowlike5, rowlike6]:
            assert collike.shape == (size, 1)
            assert collike.ndim == 2

        with pytest.raises(ValueError, match="Invalid shape"):
            data.reshape(-1, -1)
        with pytest.raises(ValueError, match="Product of shape"):
            data.reshape(len(data), 2)
