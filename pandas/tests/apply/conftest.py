import numpy as np
import pytest

from pandas import (
    DataFrame,
    Series,
)
from pandas.api.executors import BaseExecutionEngine


class MockExecutionEngine(BaseExecutionEngine):
    """
    Execution Engine to test if the execution engine interface receives and
    uses all parameters provided by the user.

    Making this engine work as the default Python engine by calling it, no extra
    functionality is implemented here.

    When testing, this will be called when this engine is provided, and then the
    same pandas.map and pandas.apply function will be called, but without engine,
    executing the default behavior from the python engine.
    """

    def map(data, func, args, kwargs, decorator, skip_na):
        kwargs_to_pass = kwargs if isinstance(data, DataFrame) else {}
        return data.map(func, na_action="ignore" if skip_na else None, **kwargs_to_pass)

    def apply(data, func, args, kwargs, decorator, axis):
        if isinstance(data, Series):
            return data.apply(func, convert_dtype=True, args=args, by_row=False)
        elif isinstance(data, DataFrame):
            return data.apply(
                func,
                axis=axis,
                raw=False,
                result_type=None,
                args=args,
                by_row="compat",
                **kwargs,
            )
        else:
            assert isinstance(data, np.ndarray)

            def wrap_function(func):
                # https://github.com/numpy/numpy/issues/8352
                def wrapper(*args, **kwargs):
                    result = func(*args, **kwargs)
                    if isinstance(result, str):
                        result = np.array(result, dtype=object)
                    return result

                return wrapper

            return np.apply_along_axis(wrap_function(func), axis, data, *args, **kwargs)


class MockEngineDecorator:
    __pandas_udf__ = MockExecutionEngine


@pytest.fixture(params=[None, MockEngineDecorator])
def engine(request):
    return request.param
