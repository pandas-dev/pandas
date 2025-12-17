import gzip
import io
import os
import subprocess
import sys
import tarfile
import textwrap
import zipfile

import numpy as np
import pytest

from pandas.compat import is_platform_windows

import pandas as pd
import pandas._testing as tm

import pandas.io.common as icom


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
def test_compression_size(obj, method, compression_only, temp_file):
    if compression_only == "tar":
        compression_only = {"method": "tar", "mode": "w:gz"}

    path = temp_file
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
def test_compression_size_fh(obj, method, compression_only, temp_file):
    path = temp_file
    with icom.get_handle(
        path,
        "w:gz" if compression_only == "tar" else "w",
        compression=compression_only,
    ) as handles:
        getattr(obj, method)(handles.handle)
        assert not handles.handle.closed
    compressed_size = os.path.getsize(path)

    # Create a new temporary file for uncompressed comparison
    path2 = temp_file.parent / f"{temp_file.stem}_uncompressed{temp_file.suffix}"
    path2.touch()
    with icom.get_handle(path2, "w", compression=None) as handles:
        getattr(obj, method)(handles.handle)
        assert not handles.handle.closed
    uncompressed_size = os.path.getsize(path2)
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
    write_method,
    write_kwargs,
    read_method,
    compression_only,
    compression_to_extension,
    temp_file,
):
    # GH22004
    input = pd.DataFrame([[1.0, 0, -4], [3.4, 5, 2]], columns=["X", "Y", "Z"])
    extension = compression_to_extension[compression_only]
    path = temp_file.parent / f"compressed{extension}"
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
    write_method,
    write_kwargs,
    read_method,
    read_kwargs,
    compression_only,
    compression_to_extension,
    temp_file,
):
    # GH22004
    input = pd.Series([0, 5, -2, 10], name="X")
    extension = compression_to_extension[compression_only]
    path = temp_file.parent / f"compressed{extension}"
    getattr(input, write_method)(path, **write_kwargs)
    if "squeeze" in read_kwargs:
        kwargs = read_kwargs.copy()
        del kwargs["squeeze"]
        output = read_method(path, compression=compression_only, **kwargs).squeeze(
            "columns"
        )
    else:
        output = read_method(path, compression=compression_only, **read_kwargs)
    tm.assert_series_equal(output, input, check_names=False)


def test_compression_warning(compression_only, temp_file):
    # Assert that passing a file object to to_csv while explicitly specifying a
    # compression protocol triggers a RuntimeWarning, as per GH21227.
    df = pd.DataFrame(
        100 * [[0.123456, 0.234567, 0.567567], [12.32112, 123123.2, 321321.2]],
        columns=["X", "Y", "Z"],
    )
    path = temp_file
    with icom.get_handle(path, "w", compression=compression_only) as handles:
        with tm.assert_produces_warning(RuntimeWarning, match="has no effect"):
            df.to_csv(handles.handle, compression=compression_only)


def test_compression_binary(compression_only, temp_file):
    """
    Binary file handles support compression.

    GH22555
    """
    df = pd.DataFrame(
        1.1 * np.arange(120).reshape((30, 4)),
        columns=pd.Index(list("ABCD")),
        index=pd.Index([f"i-{i}" for i in range(30)]),
    )

    # with a file
    path = temp_file
    with open(path, mode="wb") as file:
        df.to_csv(file, mode="wb", compression=compression_only)
        file.seek(0)  # file shouldn't be closed
    tm.assert_frame_equal(
        df, pd.read_csv(path, index_col=0, compression=compression_only)
    )

    # with BytesIO
    file = io.BytesIO()
    df.to_csv(file, mode="wb", compression=compression_only)
    file.seek(0)  # file shouldn't be closed
    tm.assert_frame_equal(
        df, pd.read_csv(file, index_col=0, compression=compression_only)
    )


def test_gzip_reproducibility_file_name(temp_file):
    """
    Gzip should create reproducible archives with mtime.

    Note: Archives created with different filenames will still be different!

    GH 28103
    """
    df = pd.DataFrame(
        1.1 * np.arange(120).reshape((30, 4)),
        columns=pd.Index(list("ABCD")),
        index=pd.Index([f"i-{i}" for i in range(30)]),
    )
    compression_options = {"method": "gzip", "mtime": 1}

    # test for filename
    path = temp_file
    df.to_csv(path, compression=compression_options)
    output = path.read_bytes()
    df.to_csv(path, compression=compression_options)
    assert output == path.read_bytes()


def test_gzip_reproducibility_file_object():
    """
    Gzip should create reproducible archives with mtime.

    GH 28103
    """
    df = pd.DataFrame(
        1.1 * np.arange(120).reshape((30, 4)),
        columns=pd.Index(list("ABCD")),
        index=pd.Index([f"i-{i}" for i in range(30)]),
    )
    compression_options = {"method": "gzip", "mtime": 1}

    # test for file object
    buffer = io.BytesIO()
    df.to_csv(buffer, compression=compression_options, mode="wb")
    output = buffer.getvalue()
    buffer = io.BytesIO()
    df.to_csv(buffer, compression=compression_options, mode="wb")
    assert output == buffer.getvalue()


@pytest.mark.single_cpu
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
    subprocess.check_output([sys.executable, "-c", code], stderr=subprocess.PIPE)


@pytest.mark.single_cpu
def test_with_missing_lzma_runtime():
    """Tests if ModuleNotFoundError is hit when calling lzma without
    having the module available.
    """
    code = textwrap.dedent(
        """
        import sys
        import pytest
        sys.modules['lzma'] = None
        import pandas as pd
        df = pd.DataFrame()
        with pytest.raises(ModuleNotFoundError, match='import of lzma'):
            df.to_csv('foo.csv', compression='xz')
        """
    )
    subprocess.check_output([sys.executable, "-c", code], stderr=subprocess.PIPE)


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
def test_gzip_compression_level(obj, method, temp_file):
    # GH33196
    path = temp_file
    getattr(obj, method)(path, compression="gzip")
    compressed_size_default = os.path.getsize(path)
    getattr(obj, method)(path, compression={"method": "gzip", "compresslevel": 1})
    compressed_size_fast = os.path.getsize(path)
    assert compressed_size_default < compressed_size_fast


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
def test_xz_compression_level_read(obj, method, temp_file):
    path = temp_file
    getattr(obj, method)(path, compression="xz")
    compressed_size_default = os.path.getsize(path)
    getattr(obj, method)(path, compression={"method": "xz", "preset": 1})
    compressed_size_fast = os.path.getsize(path)
    assert compressed_size_default < compressed_size_fast
    if method == "to_csv":
        pd.read_csv(path, compression="xz")


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
def test_bzip_compression_level(obj, method, temp_file):
    """GH33196 bzip needs file size > 100k to show a size difference between
    compression levels, so here we just check if the call works when
    compression is passed as a dict.
    """
    path = temp_file
    getattr(obj, method)(path, compression={"method": "bz2", "compresslevel": 1})


@pytest.mark.parametrize(
    "suffix,archive",
    [
        (".zip", zipfile.ZipFile),
        (".tar", tarfile.TarFile),
    ],
)
def test_empty_archive_zip(suffix, archive, temp_file):
    path = temp_file.parent / f"archive{suffix}"
    with archive(path, "w"):
        pass
    with pytest.raises(ValueError, match="Zero files found"):
        pd.read_csv(path)


def test_ambiguous_archive_zip(temp_file):
    path = temp_file.parent / "archive.zip"
    with zipfile.ZipFile(path, "w") as file:
        file.writestr("a.csv", "foo,bar")
        file.writestr("b.csv", "foo,bar")
    with pytest.raises(ValueError, match="Multiple files found in ZIP file"):
        pd.read_csv(path)


def test_ambiguous_archive_tar(tmp_path):
    csvAPath = tmp_path / "a.csv"
    with open(csvAPath, "w", encoding="utf-8") as a:
        a.write("foo,bar\n")
    csvBPath = tmp_path / "b.csv"
    with open(csvBPath, "w", encoding="utf-8") as b:
        b.write("foo,bar\n")

    tarpath = tmp_path / "archive.tar"
    with tarfile.TarFile(tarpath, "w") as tar:
        tar.add(csvAPath, "a.csv")
        tar.add(csvBPath, "b.csv")

    with pytest.raises(ValueError, match="Multiple files found in TAR archive"):
        pd.read_csv(tarpath)


def test_tar_gz_to_different_filename(temp_file):
    file = temp_file.parent / "archive.foo"
    pd.DataFrame(
        [["1", "2"]],
        columns=["foo", "bar"],
    ).to_csv(file, compression={"method": "tar", "mode": "w:gz"}, index=False)
    with gzip.open(file) as uncompressed:
        with tarfile.TarFile(fileobj=uncompressed) as archive:
            members = archive.getmembers()
            assert len(members) == 1
            content = archive.extractfile(members[0]).read().decode("utf8")

            if is_platform_windows():
                expected = "foo,bar\r\n1,2\r\n"
            else:
                expected = "foo,bar\n1,2\n"

            assert content == expected


def test_tar_no_error_on_close():
    with io.BytesIO() as buffer:
        with icom._BytesTarFile(fileobj=buffer, mode="w"):
            pass
