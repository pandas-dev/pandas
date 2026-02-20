# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import urllib.parse

import pytest
from jupyter_server.utils import url_path_join
from jupyterlab_server import LabConfig
from tornado.escape import url_escape
from traitlets import Unicode

from jupyterlab.labapp import LabApp


def mkdir(tmp_path, *parts):
    path = tmp_path.joinpath(*parts)
    if not path.exists():
        path.mkdir(parents=True)
    return path


app_settings_dir = pytest.fixture(lambda tmp_path: mkdir(tmp_path, "app_settings"))
user_settings_dir = pytest.fixture(lambda tmp_path: mkdir(tmp_path, "user_settings"))
schemas_dir = pytest.fixture(lambda tmp_path: mkdir(tmp_path, "schemas"))
workspaces_dir = pytest.fixture(lambda tmp_path: mkdir(tmp_path, "workspaces"))


@pytest.fixture
def make_lab_app(
    jp_root_dir, jp_template_dir, app_settings_dir, user_settings_dir, schemas_dir, workspaces_dir
):
    def _make_lab_app(**kwargs):
        class TestLabApp(LabApp):
            base_url = "/lab"
            extension_url = "/lab"
            default_url = Unicode("/", help="The default URL to redirect to from `/`")
            lab_config = LabConfig(
                app_name="JupyterLab Test App",
                static_dir=str(jp_root_dir),
                templates_dir=str(jp_template_dir),
                app_url="/lab",
                app_settings_dir=str(app_settings_dir),
                user_settings_dir=str(user_settings_dir),
                schemas_dir=str(schemas_dir),
                workspaces_dir=str(workspaces_dir),
            )

        app = TestLabApp()
        return app

    # Create the index files.
    index = jp_template_dir.joinpath("index.html")
    index.write_text(
        """
<!DOCTYPE html>
<html>
<head>
  <title>{{page_config['appName'] | e}}</title>
</head>
<body>
    {# Copy so we do not modify the page_config with updates. #}
    {% set page_config_full = page_config.copy() %}

    {# Set a dummy variable - we just want the side effect of the update. #}
    {% set _ = page_config_full.update(baseUrl=base_url, wsUrl=ws_url) %}

      <script id="jupyter-config-data" type="application/json">
        {{ page_config_full | tojson }}
      </script>
  <script src="{{page_config['fullStaticUrl'] | e}}/bundle.js" main="index"></script>

  <script type="text/javascript">
    /* Remove token from URL. */
    (function () {
      var parsedUrl = new URL(window.location.href);
      if (parsedUrl.searchParams.get('token')) {
        parsedUrl.searchParams.delete('token');
        window.history.replaceState({ }, '', parsedUrl.href);
      }
    })();
  </script>
</body>
</html>
"""
    )

    return _make_lab_app


@pytest.fixture
def labapp(jp_serverapp, make_lab_app):
    app = make_lab_app()
    app._link_jupyter_server_extension(jp_serverapp)
    app.initialize()
    return app


@pytest.fixture
def fetch_long(http_server_client, jp_auth_header, jp_base_url):
    """fetch fixture that handles auth, base_url, and path"""

    def client_fetch(*parts, headers=None, params=None, **kwargs):
        # Handle URL strings
        path_url = url_escape(url_path_join(*parts), plus=False)
        path_url = url_path_join(jp_base_url, path_url)
        params_url = urllib.parse.urlencode(params or {})
        url = path_url + "?" + params_url
        # Add auth keys to header
        headers = headers or {}
        headers.update(jp_auth_header)
        # Make request.
        return http_server_client.fetch(url, headers=headers, request_timeout=250, **kwargs)

    return client_fetch
