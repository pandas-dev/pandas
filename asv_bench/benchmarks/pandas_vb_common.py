from importlib import import_module
import os

import numpy as np

import pandas as pd

# Compatibility import for lib
for imp in ["pandas._libs.lib", "pandas.lib"]:
    try:
        lib = import_module(imp)
        break
    except (ImportError, TypeError, ValueError):
        pass

# Compatibility import for the testing module
try:
    import pandas._testing as tm
except ImportError:
    import pandas.util.testing as tm  # noqa: F401


numeric_dtypes = [
    np.int64,
    np.int32,
    np.uint32,
    np.uint64,
    np.float32,
    np.float64,
    np.int16,
    np.int8,
    np.uint16,
    np.uint8,
]
datetime_dtypes = [np.datetime64, np.timedelta64]
string_dtypes = [object]
try:
    extension_dtypes = [
        pd.Int8Dtype,
        pd.Int16Dtype,
        pd.Int32Dtype,
        pd.Int64Dtype,
        pd.UInt8Dtype,
        pd.UInt16Dtype,
        pd.UInt32Dtype,
        pd.UInt64Dtype,
        pd.CategoricalDtype,
        pd.IntervalDtype,
        pd.DatetimeTZDtype("ns", "UTC"),
        pd.PeriodDtype("D"),
    ]
except AttributeError:
    extension_dtypes = []


def setup(*args, **kwargs):
    # This function just needs to be imported into each benchmark file to
    # set up the random seed before each function.
    # https://asv.readthedocs.io/en/latest/writing_benchmarks.html
    np.random.seed(1234)


class BaseIO:
    """
    Base class for IO benchmarks
    """

    fname = None

    def remove(self, f):
        """Remove created files"""
        try:
            os.remove(f)
        except OSError:
            # On Windows, attempting to remove a file that is in use
            # causes an exception to be raised
            pass

    def teardown(self, *args, **kwargs):
        self.remove(self.fname)


def cow_wrapper(func):
    def new_func(*args):
        # ASV will pass us our arguments positionally
        # In the decorator, we will insert the cow parameter at the end
        # so the last arg will tell us whether CoW is set or not on the benchmark
        with pd.option_context("mode.copy_on_write", args[-1]):
            func(*args[:-1])

    return new_func


class CowDecorator:
    """
    Decorator that adds Copy on Write as a parameter to the benchmark
    """

    def __call__(self, bench_cls):
        def make_setup_wrapper(setup_f):
            def setup_wrapper(*args):
                # CoW not needed for setup, so just pass through the other parameters
                setup_f(*args[:-1])

            return setup_wrapper

        # Add CoW to the parameters
        params = getattr(bench_cls, "params", [])
        params.append([True, False])
        bench_cls.params = params
        if hasattr(bench_cls, "param_names"):
            bench_cls.param_names.append("copy_on_write")
        bench_cls.setup = make_setup_wrapper(bench_cls.setup)
        for meth in dir(bench_cls):
            if meth.startswith(("time_", "peakmem_")):
                wrapped_f = cow_wrapper(getattr(bench_cls, meth))
                setattr(bench_cls, meth, wrapped_f)
        return bench_cls
