# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import json

from jupyterlab_server.config import get_static_page_config

from jupyterlab.commands import get_app_info, lock_extension, unlock_extension
from jupyterlab.extensions.manager import PluginManager
from jupyterlab.handlers.plugin_manager_handler import PluginHandler, plugins_handler_path


def plugin_handler_labapp(jp_serverapp, make_labserver_extension_app):
    app = make_labserver_extension_app()
    app._link_jupyter_server_extension(jp_serverapp)

    # Simulate jupyter_server page config data.
    page_config = jp_serverapp.web_app.settings.setdefault("page_config_data", {})
    page_config.update(get_static_page_config(level="sys_prefix"))
    lock_rules = frozenset(
        {rule for rule, value in page_config.get("lockedExtensions", {}).items() if value}
    )

    app.handlers.extend(
        [
            (
                plugins_handler_path,
                PluginHandler,
                {
                    "manager": PluginManager(
                        ext_options={
                            "lock_rules": lock_rules,
                            "all_locked": False,
                        }
                    )
                },
            ),
        ]
    )
    return app


async def test_pluginHandler_lock_extension(jp_serverapp, jp_fetch, make_labserver_extension_app):
    extension1 = "@jupyterlab/application-extension:status"
    extension2 = "@jupyterlab/theme-dark-extension:plugin"
    info = get_app_info()
    assert info["locked"].get(extension1, False) is False
    lock_extension(extension1)
    lock_extension(extension2)
    info = get_app_info()
    assert info["locked"].get(extension1, False) is True
    assert info["locked"].get(extension2, False) is True
    # Initialize the labserver
    labapp = plugin_handler_labapp(
        jp_serverapp=jp_serverapp, make_labserver_extension_app=make_labserver_extension_app
    )
    labapp.initialize()
    # Hit and verify the pluginHandler
    response = await jp_fetch("lab", "api", "plugins", method="GET")
    payload = json.loads(response.body)
    assert response.code == 200
    assert sorted(payload["lockRules"]) == sorted([extension1, extension2])


async def test_pluginHandler_unlock_extension(jp_serverapp, jp_fetch, make_labserver_extension_app):
    extension = "@jupyterlab/application-extension:status"
    lock_extension(extension)
    info = get_app_info()
    assert info["locked"].get(extension, False) is True
    unlock_extension(extension)
    info = get_app_info()
    assert info["locked"].get(extension, False) is False
    # Initialize the labserver
    labapp = plugin_handler_labapp(
        jp_serverapp=jp_serverapp, make_labserver_extension_app=make_labserver_extension_app
    )
    labapp.initialize()
    # Hit and verify the pluginHandler
    response = await jp_fetch("lab", "api", "plugins", method="GET")
    payload = json.loads(response.body)
    assert response.code == 200
    assert payload["lockRules"] == []
