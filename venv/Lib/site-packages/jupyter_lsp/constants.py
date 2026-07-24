""" special constants used throughout jupyter_lsp
"""

# the current `entry_point` to use for python-based spec finders
EP_SPEC_V1 = "jupyter_lsp_spec_v1"

# the current `entry_point`s to use for python-based listeners
EP_LISTENER_ALL_V1 = "jupyter_lsp_listener_all_v1"
EP_LISTENER_CLIENT_V1 = "jupyter_lsp_listener_client_v1"
EP_LISTENER_SERVER_V1 = "jupyter_lsp_listener_server_v1"

# jupyter*config.d where language_servers can be defined
APP_CONFIG_D_SECTIONS = ["_", "_notebook_", "_server_"]
