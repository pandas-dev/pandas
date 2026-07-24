from typing import Final

from django.apps import AppConfig

class ChannelsConfig(AppConfig):
    name: Final = "channels"
    verbose_name: str = "Channels"
