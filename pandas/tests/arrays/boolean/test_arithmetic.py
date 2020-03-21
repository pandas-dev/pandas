import numpy as np
import pytest

import pandas as pd
from pandas.tests.extension.base import BaseOpsUtil


@pytest.fixture
def data():
    return pd.array(
        [True, False] * 4 + [np.nan] + [True, False] * 44 + [np.nan] + [True, False],
        dtype="boolean",
    )


class TestArithmeticOps(BaseOpsUtil):
    def test_error(self, data, all_arithmetic_operators):
        # invalid ops

        op = all_arithmetic_operators
        s = pd.Series(data)
        ops = getattr(s, op)
        opa = getattr(data, op)

        # invalid scalars
        with pytest.raises(TypeError):
            ops("foo")
        with pytest.raises(TypeError):
            ops(pd.Timestamp("20180101"))

        # invalid array-likes
        if op not in ("__mul__", "__rmul__"):
            # TODO(extension) numpy's mul with object array sees booleans as numbers
            with pytest.raises(TypeError):
                ops(pd.Series("foo", index=s.index))

        # 2d
        result = opa(pd.DataFrame({"A": s}))
        assert result is NotImplemented

        with pytest.raises(NotImplementedError):
            opa(np.arange(len(s)).reshape(-1, len(s)))
