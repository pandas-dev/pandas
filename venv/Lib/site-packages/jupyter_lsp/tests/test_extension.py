import os

import pytest


def test_serverextension_path(app):
    import jupyter_lsp

    paths = jupyter_lsp._jupyter_server_extension_paths()
    for path in paths:
        assert __import__(path["module"])


def test_serverextension(app):
    app.initialize(
        ["--ServerApp.jpserver_extensions={'jupyter_lsp.serverextension': True}"]
    )
    assert app.language_server_manager
    found_lsp = False
    for r in app.web_app.default_router.rules:
        for rr in r.target.rules:
            if "/lsp/" in str(rr.matcher.regex):
                found_lsp = True

    assert found_lsp, "apparently didn't install the /lsp/ route"


def test_default_virtual_documents_dir(app):
    app.initialize(
        ["--ServerApp.jpserver_extensions={'jupyter_lsp.serverextension': True}"]
    )
    assert app.language_server_manager.virtual_documents_dir == ".virtual_documents"


def test_virtual_documents_dir_config(app):
    custom_dir = ".custom_virtual_dir"
    app.initialize(
        [
            "--ServerApp.jpserver_extensions={'jupyter_lsp.serverextension': True}",
            "--ServerApp.LanguageServerManager.virtual_documents_dir=" + custom_dir,
        ]
    )
    assert app.language_server_manager.virtual_documents_dir == custom_dir


def test_virtual_documents_dir_env(app):
    os.environ["JP_LSP_VIRTUAL_DIR"] = custom_dir = ".custom_virtual_dir"
    app.initialize(
        ["--ServerApp.jpserver_extensions={'jupyter_lsp.serverextension': True}"]
    )
    assert app.language_server_manager.virtual_documents_dir == custom_dir
    page_config = app.web_app.settings["page_config_data"]
    assert page_config["virtualDocumentsUri"].endswith(custom_dir)


@pytest.mark.parametrize("invalid_path", ["", "."])
def test_virtual_documents_dir_env_empty(app, invalid_path):
    os.environ["JP_LSP_VIRTUAL_DIR"] = invalid_path
    app.initialize(
        ["--ServerApp.jpserver_extensions={'jupyter_lsp.serverextension': True}"]
    )
    page_config = app.web_app.settings["page_config_data"]
    assert page_config["virtualDocumentsUri"].endswith("/.virtual_documents")
