import contextlib
import os
import subprocess
import textwrap
import warnings

import pytest

import pandas as pd
import pandas.util.testing as tm

import pandas.io.common as icom


@contextlib.contextmanager
def catch_to_csv_depr():
    # Catching warnings because Series.to_csv has
    # been deprecated. Remove this context when
    # Series.to_csv has been aligned.

    with warnings.catch_warnings(record=True):
        warnings.simplefilter("ignore", FutureWarning)
        yield


@pytest.mark.parametrize(
    "obj",
    [
        pd.DataFrame(
            100 * [[0.123456, 0.234567, 0.567567], [12.32112, 123123.2, 321321.2]],
            columns=["X", "Y", "Z"],
        ),
        pd.Series(100 * [0.123456, 0.234567, 0.567567], name="X"),
    ],
)
@pytest.mark.parametrize("method", ["to_pickle", "to_json", "to_csv"])
def test_compression_size(obj, method, compression_only):
    with tm.ensure_clean() as path:
        with catch_to_csv_depr():
            getattr(obj, method)(path, compression=compression_only)
            compressed_size = os.path.getsize(path)
            getattr(obj, method)(path, compression=None)
            uncompressed_size = os.path.getsize(path)
            assert uncompressed_size > compressed_size


@pytest.mark.parametrize(
    "obj",
    [
        pd.DataFrame(
            100 * [[0.123456, 0.234567, 0.567567], [12.32112, 123123.2, 321321.2]],
            columns=["X", "Y", "Z"],
        ),
        pd.Series(100 * [0.123456, 0.234567, 0.567567], name="X"),
    ],
)
@pytest.mark.parametrize("method", ["to_csv", "to_json"])
def test_compression_size_fh(obj, method, compression_only):
    with tm.ensure_clean() as path:
        f, handles = icom._get_handle(path, "w", compression=compression_only)
        with catch_to_csv_depr():
            with f:
                getattr(obj, method)(f)
                assert not f.closed
            assert f.closed
            compressed_size = os.path.getsize(path)
    with tm.ensure_clean() as path:
        f, handles = icom._get_handle(path, "w", compression=None)
        with catch_to_csv_depr():
            with f:
                getattr(obj, method)(f)
                assert not f.closed
        assert f.closed
        uncompressed_size = os.path.getsize(path)
        assert uncompressed_size > compressed_size


@pytest.mark.parametrize(
    "write_method, write_kwargs, read_method",
    [
        ("to_csv", {"index": False}, pd.read_csv),
        ("to_json", {}, pd.read_json),
        ("to_pickle", {}, pd.read_pickle),
    ],
)
def test_dataframe_compression_defaults_to_infer(
    write_method, write_kwargs, read_method, compression_only
):
    # GH22004
    input = pd.DataFrame([[1.0, 0, -4], [3.4, 5, 2]], columns=["X", "Y", "Z"])
    extension = icom._compression_to_extension[compression_only]
    with tm.ensure_clean("compressed" + extension) as path:
        getattr(input, write_method)(path, **write_kwargs)
        output = read_method(path, compression=compression_only)
    tm.assert_frame_equal(output, input)


@pytest.mark.parametrize(
    "write_method,write_kwargs,read_method,read_kwargs",
    [
        ("to_csv", {"index": False, "header": True}, pd.read_csv, {"squeeze": True}),
        ("to_json", {}, pd.read_json, {"typ": "series"}),
        ("to_pickle", {}, pd.read_pickle, {}),
    ],
)
def test_series_compression_defaults_to_infer(
    write_method, write_kwargs, read_method, read_kwargs, compression_only
):
    # GH22004
    input = pd.Series([0, 5, -2, 10], name="X")
    extension = icom._compression_to_extension[compression_only]
    with tm.ensure_clean("compressed" + extension) as path:
        getattr(input, write_method)(path, **write_kwargs)
        output = read_method(path, compression=compression_only, **read_kwargs)
    tm.assert_series_equal(output, input, check_names=False)


def test_compression_warning(compression_only):
    # Assert that passing a file object to to_csv while explicitly specifying a
    # compression protocol triggers a RuntimeWarning, as per GH21227.
    df = pd.DataFrame(
        100 * [[0.123456, 0.234567, 0.567567], [12.32112, 123123.2, 321321.2]],
        columns=["X", "Y", "Z"],
    )
    with tm.ensure_clean() as path:
        f, handles = icom._get_handle(path, "w", compression=compression_only)
        with tm.assert_produces_warning(RuntimeWarning, check_stacklevel=False):
            with f:
                df.to_csv(f, compression=compression_only)


def test_with_missing_lzma():
    """Tests if import pandas works when lzma is not present."""
    # https://github.com/pandas-dev/pandas/issues/27575
    code = textwrap.dedent(
        """\
        import sys
        sys.modules['lzma'] = None
        import pandas
        """
    )
    subprocess.check_output(["python", "-c", code])


def test_with_missing_lzma_runtime():
    """Tests if RuntimeError is hit when calling lzma without
    having the module available."""
    code = textwrap.dedent(
        """
        import sys
        import pytest
        sys.modules['lzma'] = None
        import pandas
        df = pandas.DataFrame()
        with pytest.raises(RuntimeError, match='lzma module'):
            df.to_csv('foo.csv', compression='xz')
        """
    )
    subprocess.check_output(["python", "-c", code])
