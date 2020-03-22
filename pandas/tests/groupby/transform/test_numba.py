import numpy as np
import pytest

import pandas.util._test_decorators as td

from pandas import DataFrame, Series
import pandas._testing as tm

@td.skip_if_no("numba", "0.46.0")
def test_correct_function_signature():
    pass