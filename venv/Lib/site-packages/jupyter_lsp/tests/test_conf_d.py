def test_conf_d_language_servers(echo_conf_json, handlers, app_config_d):
    (app_config_d / "echo.json").write_text(echo_conf_json)
    handler, ws_handler = handlers
    manager = handler.manager
    manager.initialize()
    assert "_echo_" in [*manager.language_servers]
