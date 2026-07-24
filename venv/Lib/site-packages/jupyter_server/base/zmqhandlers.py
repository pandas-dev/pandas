"""This module is deprecated in Jupyter Server 2.0"""

# Raise a warning that this module is deprecated.
import warnings

from tornado.websocket import WebSocketHandler

from jupyter_server.base.websocket import WebSocketMixin
from jupyter_server.services.kernels.connection.base import (
    deserialize_binary_message,
    deserialize_msg_from_ws_v1,
    serialize_binary_message,
    serialize_msg_to_ws_v1,
)

warnings.warn(
    "jupyter_server.base.zmqhandlers module is deprecated in Jupyter Server 2.0",
    DeprecationWarning,
    stacklevel=2,
)
