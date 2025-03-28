import pathlib
import platform
import sys

import pytest

from ..paths import file_uri_to_path, is_relative, normalized_uri

WIN = platform.system() == "Windows"
HOME = pathlib.Path("~").expanduser()
PY35 = sys.version_info[:2] == (3, 5)


@pytest.mark.skipif(WIN, reason="can't test POSIX paths on Windows")
@pytest.mark.parametrize("root_dir, expected_root_uri", [["~", HOME.as_uri()]])
def test_normalize_posix_path_home(root_dir, expected_root_uri):  # pragma: no cover
    assert normalized_uri(root_dir) == expected_root_uri


@pytest.mark.skipif(WIN, reason="can't test POSIX paths on Windows")
@pytest.mark.parametrize(
    "root, path",
    [["~", "~/a"], ["~", "~/a/../b/"], ["/", "/"], ["/a", "/a/b"], ["/a", "/a/b/../c"]],
)
def test_is_relative_ok(root, path):
    assert is_relative(root, path)


@pytest.mark.skipif(WIN, reason="can't test POSIX paths on Windows")
@pytest.mark.parametrize(
    "root, path",
    [
        ["~", "~/.."],
        ["~", "/"],
        ["/a", "/"],
        ["/a/b", "/a"],
        ["/a/b", "/a/b/.."],
        ["/a", "/a/../b"],
        ["/a", "a//"],
    ],
)
def test_is_relative_not_ok(root, path):
    assert not is_relative(root, path)


@pytest.mark.skipif(not WIN, reason="can't test Windows paths on POSIX")
@pytest.mark.parametrize(
    "root, path",
    [
        ["c:\\Users\\user1", "c:\\Users\\"],
        ["c:\\Users\\user1", "d:\\"],
        ["c:\\Users", "c:\\Users\\.."],
    ],
)
def test_is_relative_not_ok_win(root, path):
    assert not is_relative(root, path)


@pytest.mark.skipif(PY35, reason="can't test non-existent paths on py35")
@pytest.mark.skipif(WIN, reason="can't test POSIX paths on Windows")
@pytest.mark.parametrize(
    "root_dir, expected_root_uri",
    [
        # probably need to try some other things
        [str(HOME / "foo"), (HOME / "foo").as_uri()]
    ],
)
def test_normalize_posix_path_home_subdir(
    root_dir, expected_root_uri
):  # pragma: no cover
    assert normalized_uri(root_dir) == expected_root_uri


@pytest.mark.skipif(not WIN, reason="can't test Windows paths on POSIX")
@pytest.mark.parametrize(
    "root_dir, expected_root_uri",
    [
        ["c:\\Users\\user1", "file:///c:/Users/user1"],
        ["C:\\Users\\user1", "file:///c:/Users/user1"],
        ["//VBOXSVR/shared-folder", "file://vboxsvr/shared-folder/"],
    ],
)
def test_normalize_windows_path_case(root_dir, expected_root_uri):  # pragma: no cover
    try:
        normalized = normalized_uri(root_dir)
    except FileNotFoundError as err:
        if sys.version_info >= (3, 10):
            # apparently, this triggers resolving the path on win/py3.10
            return
        raise err

    assert normalized == expected_root_uri


@pytest.mark.skipif(WIN, reason="can't test POSIX paths on Windows")
@pytest.mark.parametrize(
    "file_uri, expected_posix_path",
    [
        ["file:///C:/Windows/System32/Drivers/etc", "/C:/Windows/System32/Drivers/etc"],
        ["file:///C:/some%20dir/some%20file.txt", "/C:/some dir/some file.txt"],
        ["file:///home/user/some%20file.txt", "/home/user/some file.txt"],
    ],
)
def test_file_uri_to_path_posix(file_uri, expected_posix_path):  # pragma: no cover
    assert file_uri_to_path(file_uri) == expected_posix_path


@pytest.mark.skipif(not WIN, reason="can't test Windows paths on POSIX")
@pytest.mark.parametrize(
    "file_uri, expected_windows_path",
    [
        # https://github.com/jupyter-lsp/jupyterlab-lsp/pull/305#issuecomment-665996145
        ["file:///C:/Windows/System32/Drivers/etc", "C:/Windows/System32/Drivers/etc"],
        ["file:///C:/some%20dir/some%20file.txt", "C:/some dir/some file.txt"],
    ],
)
def test_file_uri_to_path_windows(file_uri, expected_windows_path):  # pragma: no cover
    assert file_uri_to_path(file_uri) == expected_windows_path
