"""ApiGatewayManagementApiBackend class with methods for supported APIs."""

from collections import defaultdict
from typing import Any, Dict

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.utils import unix_time


class Connection:
    def __init__(self) -> None:
        self.connected_at = unix_time()
        self.source_ip = "192.168.0.1"
        self.user_agent = "Moto Mocks"
        self.data = b""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "connectedAt": self.connected_at,
            "lastActiveAt": unix_time(),
            "identity": {
                "sourceIp": self.source_ip,
                "userAgent": self.user_agent,
            },
        }


class ApiGatewayManagementApiBackend(BaseBackend):
    """
    Connecting to this API in ServerMode/Docker requires Python >= 3.8 and an up-to-date `werkzeug` version (>=2.3.x)
    """

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.connections: Dict[str, Connection] = defaultdict(Connection)

    def delete_connection(self, connection_id: str) -> None:
        self.connections.pop(connection_id, None)

    def get_connection(self, connection_id: str) -> Connection:
        return self.connections[connection_id]

    def post_to_connection(self, data: bytes, connection_id: str) -> None:
        cnctn = self.get_connection(connection_id)
        cnctn.data += data


apigatewaymanagementapi_backends = BackendDict(
    ApiGatewayManagementApiBackend, "apigateway"
)
