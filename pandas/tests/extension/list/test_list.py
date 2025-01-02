import pytest

import pandas as pd
from pandas.core.arrays.list_ import (
    ListArray,
    ListDtype,
)
from pandas.tests.extension.base.accumulate import BaseAccumulateTests
from pandas.tests.extension.base.casting import BaseCastingTests
from pandas.tests.extension.base.constructors import BaseConstructorsTests
from pandas.tests.extension.base.dim2 import (  # noqa: F401
    Dim2CompatTests,
    NDArrayBacked2DTests,
)
from pandas.tests.extension.base.dtype import BaseDtypeTests
from pandas.tests.extension.base.getitem import BaseGetitemTests
from pandas.tests.extension.base.groupby import BaseGroupbyTests
from pandas.tests.extension.base.index import BaseIndexTests
from pandas.tests.extension.base.interface import BaseInterfaceTests
from pandas.tests.extension.base.io import BaseParsingTests
from pandas.tests.extension.base.methods import BaseMethodsTests
from pandas.tests.extension.base.missing import BaseMissingTests
from pandas.tests.extension.base.ops import (  # noqa: F401
    BaseArithmeticOpsTests,
    BaseComparisonOpsTests,
    BaseOpsUtil,
    BaseUnaryOpsTests,
)
from pandas.tests.extension.base.printing import BasePrintingTests
from pandas.tests.extension.base.reduce import BaseReduceTests
from pandas.tests.extension.base.reshaping import BaseReshapingTests
from pandas.tests.extension.base.setitem import BaseSetitemTests



@pytest.fixture
def dtype():
    return ListDtype()


@pytest.fixture
def data():
    """Length-100 ListArray for semantics test."""
    # TODO: make better random data
    data = [list("a"), list("ab"), list("abc")] * 33 + [None]
    return ListArray(data)


class TestListArray(
    BaseAccumulateTests,
    #BaseCastingTests,
    BaseConstructorsTests,
    #BaseDtypeTests,
    #BaseGetitemTests,
    #BaseGroupbyTests,
    BaseIndexTests,
    #BaseInterfaceTests,
    BaseParsingTests,
    #BaseMethodsTests,
    #BaseMissingTests,
    #BaseArithmeticOpsTests,
    #BaseComparisonOpsTests,
    #BaseUnaryOpsTests,
    #BasePrintingTests,
    BaseReduceTests,
    #BaseReshapingTests,
    #BaseSetitemTests,
    Dim2CompatTests,
):
    ...


def test_to_csv(data):
    # https://github.com/pandas-dev/pandas/issues/28840
    # array with list-likes fail when doing astype(str) on the numpy array
    # which was done in get_values_for_csv
    df = pd.DataFrame({"a": data})
    res = df.to_csv()
    assert str(data[0]) in res
