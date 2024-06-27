from abc import ABC, abstractmethod
from typing import Any, List


class KernelWebsocketConnectionABC(ABC):
    """
    This class defines a minimal interface that should
    be used to bridge the connection between Jupyter
    Server's websocket API and a kernel's ZMQ socket
    interface.
    """

    websocket_handler: Any

    @abstractmethod
    async def connect(self):
        """Connect the kernel websocket to the kernel ZMQ connections"""

    @abstractmethod
    async def disconnect(self):
        """Disconnect the kernel websocket from the kernel ZMQ connections"""

    @abstractmethod
    def handle_incoming_message(self, incoming_msg: str) -> None:
        """Broker the incoming websocket message to the appropriate ZMQ channel."""

    @abstractmethod
    def handle_outgoing_message(self, stream: str, outgoing_msg: List[Any]) -> None:
        """Broker outgoing ZMQ messages to the kernel websocket."""
