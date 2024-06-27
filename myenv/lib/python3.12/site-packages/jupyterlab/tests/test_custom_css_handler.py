# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import os

import pytest

CUSTOM_CSS = """body #top-panel-wrapper,
#jp-top-bar {
  background-color: #aecad4 !important;
}

body h1 {
  font-size: 22px;
  margin-bottom: 40px;
  color: #10929e;
  text-decoration: underline;
}"""


@pytest.fixture
def jp_server_config(jp_server_config, tmp_path):
    config = jp_server_config.copy()
    config["LabApp"]["custom_css"] = True
    return config


async def test_CustomCssHandler(tmp_path, jp_serverapp, labserverapp, jp_fetch):
    custom_path = tmp_path / "config" / "custom"
    # Check we are placing the custom.css file in the appropriate folder
    assert str(custom_path) in jp_serverapp.web_app.settings["static_custom_path"]
    custom_path.mkdir(parents=True, exist_ok=True)
    (custom_path / "custom.css").write_text(CUSTOM_CSS)
    response = await jp_fetch("custom", "custom.css", method="GET")

    assert response.code == 200
    assert response.body.decode().replace(os.linesep, "\n") == CUSTOM_CSS
