import collections
from functools import partial
import string
import subprocess
import sys

import numpy as np
import pytest

from pandas.compat import WASM

import pandas as pd
from pandas import Series
import pandas._testing as tm
from pandas.core import ops
import pandas.core.common as com
from pandas.util.version import Version


class TestGetCallableName:
    def fn(self, x):
        return x

    partial1 = partial(fn)
    partial2 = partial(partial1)
    lambda_ = lambda x: x

    class SomeCall:
        def __call__(self):
            # This shouldn't actually get called below; SomeCall.__init__
            #  should.
            raise NotImplementedError

    @pytest.mark.parametrize(
        "func, expected",
        [
            (fn, "fn"),
            (partial1, "fn"),
            (partial2, "fn"),
            (lambda_, "<lambda>"),
            (SomeCall(), "SomeCall"),
            (1, None),
        ],
    )
    def test_get_callable_name(self, func, expected):
        assert com.get_callable_name(func) == expected


class TestRandomState:
    def test_seed(self):
        seed = 5
        assert com.random_state(seed).uniform() == np.random.RandomState(seed).uniform()

    def test_object(self):
        seed = 10
        state_obj = np.random.RandomState(seed)
        assert (
            com.random_state(state_obj).uniform()
            == np.random.RandomState(seed).uniform()
        )

    def test_default(self):
        assert com.random_state() is np.random

    def test_array_like(self):
        state = np.random.default_rng(None).integers(0, 2**31, size=624, dtype="uint32")
        assert (
            com.random_state(state).uniform() == np.random.RandomState(state).uniform()
        )

    def test_bit_generators(self):
        seed = 3
        assert (
            com.random_state(np.random.MT19937(seed)).uniform()
            == np.random.RandomState(np.random.MT19937(seed)).uniform()
        )

        seed = 11
        assert (
            com.random_state(np.random.PCG64(seed)).uniform()
            == np.random.RandomState(np.random.PCG64(seed)).uniform()
        )

    @pytest.mark.parametrize("state", ["test", 5.5])
    def test_error(self, state):
        msg = (
            "random_state must be an integer, array-like, a BitGenerator, Generator, "
            "a numpy RandomState, or None"
        )
        with pytest.raises(ValueError, match=msg):
            com.random_state(state)


@pytest.mark.parametrize("args, expected", [((1, 2, None), True), ((1, 2, 3), False)])
def test_any_none(args, expected):
    assert com.any_none(*args) is expected


@pytest.mark.parametrize(
    "args, expected",
    [((1, 2, 3), True), ((1, 2, None), False), ((None, None, None), False)],
)
def test_all_not_none(args, expected):
    assert com.all_not_none(*args) is expected


@pytest.mark.parametrize(
    "left, right, expected",
    [
        (Series([1], name="x"), Series([2], name="x"), "x"),
        (Series([1], name="x"), Series([2], name="y"), None),
        (Series([1]), Series([2], name="x"), None),
        (Series([1], name="x"), Series([2]), None),
        (Series([1], name="x"), [2], "x"),
        ([1], Series([2], name="y"), "y"),
        # matching NAs
        (Series([1], name=np.nan), pd.Index([], name=np.nan), np.nan),
        (Series([1], name=np.nan), pd.Index([], name=pd.NaT), None),
        (Series([1], name=pd.NA), pd.Index([], name=pd.NA), pd.NA),
        # tuple name GH#39757
        (
            Series([1], name=np.int64(1)),
            pd.Index([], name=(np.int64(1), np.int64(2))),
            None,
        ),
        (
            Series([1], name=(np.int64(1), np.int64(2))),
            pd.Index([], name=(np.int64(1), np.int64(2))),
            (np.int64(1), np.int64(2)),
        ),
        pytest.param(
            Series([1], name=(np.float64("nan"), np.int64(2))),
            pd.Index([], name=(np.float64("nan"), np.int64(2))),
            (np.float64("nan"), np.int64(2)),
            marks=pytest.mark.xfail(
                reason="Not checking for matching NAs inside tuples."
            ),
        ),
    ],
)
def test_maybe_match_name(left, right, expected):
    res = ops.common._maybe_match_name(left, right)
    assert res is expected or res == expected


@pytest.mark.parametrize(
    "into, msg",
    [
        (
            # uninitialized defaultdict
            collections.defaultdict,
            r"to_dict\(\) only accepts initialized defaultdicts",
        ),
        (
            # non-mapping subtypes,, instance
            [],
            "unsupported type: <class 'list'>",
        ),
        (
            # non-mapping subtypes, class
            list,
            "unsupported type: <class 'list'>",
        ),
    ],
)
def test_standardize_mapping_type_error(into, msg):
    with pytest.raises(TypeError, match=msg):
        com.standardize_mapping(into)


def test_standardize_mapping():
    fill = {"bad": "data"}
    assert com.standardize_mapping(fill) == dict

    # Convert instance to type
    assert com.standardize_mapping({}) == dict

    dd = collections.defaultdict(list)
    assert isinstance(com.standardize_mapping(dd), partial)


def test_git_version():
    # GH 21295
    git_version = pd.__git_version__
    assert len(git_version) == 40
    assert all(c in string.hexdigits for c in git_version)


def test_version_tag():
    version = Version(pd.__version__)
    try:
        version > Version("0.0.1")
    except TypeError as err:
        raise ValueError(
            "No git tags exist, please sync tags between upstream and your repo"
        ) from err


@pytest.mark.parametrize("obj", [obj for obj in pd.__dict__.values() if callable(obj)])
def test_serializable(obj, temp_file):
    # GH 35611
    unpickled = tm.round_trip_pickle(obj, temp_file)
    assert type(obj) == type(unpickled)


class TestIsBoolIndexer:
    def test_non_bool_array_with_na(self):
        # in particular, this should not raise
        arr = np.array(["A", "B", np.nan], dtype=object)
        assert not com.is_bool_indexer(arr)

    def test_list_subclass(self):
        # GH#42433

        class MyList(list):
            pass

        val = MyList(["a"])

        assert not com.is_bool_indexer(val)

        val = MyList([True])
        assert com.is_bool_indexer(val)

    def test_frozenlist(self):
        # GH#42461
        data = {"col1": [1, 2], "col2": [3, 4]}
        df = pd.DataFrame(data=data)

        frozen = df.index.names[1:]
        assert not com.is_bool_indexer(frozen)

        result = df[frozen]
        expected = df[[]]
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("scalar", [1, True])
    def test_numpyextensionarray(self, scalar):
        # GH 63391
        arr = pd.arrays.NumpyExtensionArray(np.array([scalar]))
        assert com.is_bool_indexer(arr) is isinstance(scalar, bool)


@pytest.mark.parametrize("with_exception", [True, False])
def test_temp_setattr(with_exception):
    # GH#45954
    ser = Series(dtype=object)
    ser.name = "first"
    # Raise a ValueError in either case to satisfy pytest.raises
    match = "Inside exception raised" if with_exception else "Outside exception raised"
    with pytest.raises(ValueError, match=match):
        with com.temp_setattr(ser, "name", "second"):
            assert ser.name == "second"
            if with_exception:
                raise ValueError("Inside exception raised")
        raise ValueError("Outside exception raised")
    assert ser.name == "first"


@pytest.mark.skipif(WASM, reason="Can't start subprocesses in WASM")
@pytest.mark.single_cpu
def test_str_size():
    # GH#21758
    a = "a"
    expected = sys.getsizeof(a)
    pyexe = sys.executable.replace("\\", "/")
    call = [
        pyexe,
        "-c",
        "a='a';import sys;sys.getsizeof(a);import pandas;print(sys.getsizeof(a));",
    ]
    result = subprocess.check_output(call).decode()[-4:-1].strip("\n")
    assert int(result) == int(expected)
