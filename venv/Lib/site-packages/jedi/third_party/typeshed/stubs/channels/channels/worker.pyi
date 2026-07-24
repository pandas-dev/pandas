from asgiref.server import StatelessServer
from channels.layers import BaseChannelLayer
from channels.utils import _ChannelApplication

class Worker(StatelessServer):
    channels: list[str]
    channel_layer: BaseChannelLayer

    def __init__(
        self, application: _ChannelApplication, channels: list[str], channel_layer: BaseChannelLayer, max_applications: int = 1000
    ) -> None: ...
    async def handle(self) -> None: ...
    async def listener(self, channel: str) -> None: ...
