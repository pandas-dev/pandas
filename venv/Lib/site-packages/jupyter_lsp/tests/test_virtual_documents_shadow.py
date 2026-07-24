import logging
from pathlib import Path
from types import SimpleNamespace
from typing import List

import pytest

from jupyter_lsp import LanguageServerManager

from ..virtual_documents_shadow import (
    EditableFile,
    ShadowFilesystemError,
    extract_or_none,
    setup_shadow_filesystem,
)


@pytest.mark.asyncio
async def test_read(tmp_path):
    path = tmp_path / "existing.py"
    path.write_text("a\ntest")

    editable_file = EditableFile(path)

    await editable_file.read()

    assert editable_file.lines == ["a", "test"]


@pytest.mark.asyncio
async def test_read_missing(tmp_path):
    path = tmp_path / "missing.py"
    missing_file = EditableFile(path)

    await missing_file.read()

    assert missing_file.lines == [""]


@pytest.mark.asyncio
async def test_apply_change(tmp_path):
    # inserting text
    path = tmp_path / "test.py"
    editable_file = EditableFile(path)
    await editable_file.read()

    editable_file.apply_change("new\ntext", **editable_file.full_range)
    assert editable_file.lines == ["new", "text"]

    # modifying a range
    editable_file.apply_change(
        "ves", start={"line": 1, "character": 0}, end={"line": 1, "character": 3}
    )
    assert editable_file.lines == ["new", "vest"]

    editable_file.apply_change("", **editable_file.full_range)
    assert editable_file.lines == [""]


def test_extract_or_none():
    obj = {"nested": {"value": 1}}
    assert extract_or_none(obj, ["nested"]) == {"value": 1}
    assert extract_or_none(obj, ["nested", "value"]) == 1
    assert extract_or_none(obj, ["missing", "value"]) is None


def did_open(uri, text):
    return {
        "method": "textDocument/didOpen",
        "params": {"textDocument": {"uri": uri, "text": text}},
    }


def did_change(uri, changes: List):
    return {
        "method": "textDocument/didChange",
        "params": {"textDocument": {"uri": uri}, "contentChanges": changes},
    }


def did_save_with_text(uri, text):
    return {
        "method": "textDocument/didSave",
        "params": {"textDocument": {"uri": uri, "text": text}},
    }


def did_save_without_text(uri):
    return {"method": "textDocument/didSave", "params": {"textDocument": {"uri": uri}}}


@pytest.fixture
def shadow_path(tmpdir):
    return str(tmpdir.mkdir(".virtual_documents"))


@pytest.fixture
def manager():
    manager = LanguageServerManager()
    manager.language_servers = {
        "python-lsp-server": {
            "requires_documents_on_disk": True,
            "argv": [],
            "languages": ["python"],
            "version": 2,
        }
    }
    return manager


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "message_func, content, expected_content",
    [
        [did_open, "content\nof\nopened\nfile", "content\nof\nopened\nfile"],
        [did_change, [{"text": "content after change"}], "content after change"],
        [did_save_with_text, "content at save", "content at save"],
    ],
)
async def test_shadow(shadow_path, message_func, content, expected_content, manager):
    shadow = setup_shadow_filesystem(Path(shadow_path).as_uri())
    ok_file_path = Path(shadow_path) / "test.py"

    message = message_func(ok_file_path.as_uri(), content)
    result = await shadow("client", message, "python-lsp-server", manager)
    assert isinstance(result, str)

    with open(ok_file_path) as f:
        assert f.read() == expected_content


@pytest.mark.asyncio
async def test_no_shadow_for_well_behaved_server(
    shadow_path,
):
    """We call server well behaved when it does not require a disk copy"""
    shadow_path_for_well = Path(shadow_path) / "well"
    shadow = setup_shadow_filesystem(Path(shadow_path_for_well).as_uri())
    ok_file_path = Path(shadow_path_for_well) / "test.py"

    manager = SimpleNamespace(
        language_servers={"python-lsp-server": {"requires_documents_on_disk": False}}
    )

    message = did_open(ok_file_path.as_uri(), "content\nof\nopened\nfile")
    result = await shadow("client", message, "python-lsp-server", manager)
    # should short-circuit for well behaved server
    assert result is None
    # should not create the directory
    assert not shadow_path_for_well.exists()


@pytest.mark.asyncio
async def test_shadow_created_for_ill_behaved_server(
    shadow_path,
):
    shadow_path_for_ill = Path(shadow_path) / "ill"
    shadow = setup_shadow_filesystem(shadow_path_for_ill.as_uri())
    ok_file_path = Path(shadow_path_for_ill) / "test.py"

    manager = SimpleNamespace(
        language_servers={"python-lsp-server": {"requires_documents_on_disk": True}}
    )

    message = did_open(ok_file_path.as_uri(), "content\nof\nopened\nfile")
    result = await shadow("client", message, "python-lsp-server", manager)
    assert result is not None
    # should create the directory at given path
    assert shadow_path_for_ill.exists()
    assert shadow_path_for_ill.is_dir()


@pytest.mark.asyncio
async def test_shadow_failures(shadow_path, manager):
    shadow = setup_shadow_filesystem(Path(shadow_path).as_uri())
    ok_file_uri = (Path(shadow_path) / "test.py").as_uri()

    def run_shadow(message):
        return shadow("client", message, "python-lsp-server", manager)

    # missing textDocument
    with pytest.raises(ShadowFilesystemError, match="Could not get textDocument from"):
        await run_shadow({"method": "textDocument/didChange"})

    # missing URI
    with pytest.raises(ShadowFilesystemError, match="Could not get URI from"):
        await run_shadow(
            {"method": "textDocument/didChange", "params": {"textDocument": {}}}
        )

    # should ignore other methods
    result = await run_shadow({"method": "textDocument/completion"})
    assert result is None

    # should NOT intercept (nor shadow) files from location other than shadow_path
    result = await run_shadow(did_open("file:///other/path.py", "content"))
    assert result is None

    # should fail silently on missing text in didSave
    result = await run_shadow(did_save_without_text(ok_file_uri))
    assert result is None

    # should raise on missing changes in didChange
    with pytest.raises(ShadowFilesystemError, match=".* is missing contentChanges"):
        await run_shadow(
            {
                "method": "textDocument/didChange",
                "params": {"textDocument": {"uri": ok_file_uri}},
            }
        )


@pytest.mark.asyncio
async def test_shadow_traversal(shadow_path, manager):
    file_beyond_shadow_root_uri = (Path(shadow_path) / ".." / "test.py").as_uri()

    shadow = setup_shadow_filesystem(Path(shadow_path).as_uri())

    def run_shadow(message):
        return shadow("client", message, "python-lsp-server", manager)

    with pytest.raises(
        ShadowFilesystemError, match="is not relative to shadow filesystem root"
    ):
        await run_shadow(did_open(file_beyond_shadow_root_uri, "content"))


@pytest.fixture
def forbidden_shadow_path(tmpdir):
    path = Path(tmpdir) / "no_permission_dir"
    path.mkdir()
    path.chmod(0o000)

    yield path

    # re-adjust the permissions, see https://github.com/pytest-dev/pytest/issues/7821
    path.chmod(0o755)


@pytest.mark.asyncio
async def test_io_failure(forbidden_shadow_path, manager, caplog):
    file_uri = (forbidden_shadow_path / "test.py").as_uri()

    shadow = setup_shadow_filesystem(forbidden_shadow_path.as_uri())

    def send_change():
        message = did_open(file_uri, "content")
        return shadow("client", message, "python-lsp-server", manager)

    with caplog.at_level(logging.WARNING):
        assert await send_change() is None
        assert await send_change() is None
    # no message should be emitted during the first two attempts
    assert caplog.text == ""

    # a warning should be emitted on third failure
    with caplog.at_level(logging.WARNING):
        assert await send_change() is None
    assert "initialization of shadow filesystem failed three times" in caplog.text
    assert "PermissionError" in caplog.text
    caplog.clear()

    # no message should be emitted in subsequent attempts
    with caplog.at_level(logging.WARNING):
        assert await send_change() is None
    assert caplog.text == ""
